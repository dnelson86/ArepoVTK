/*
 * arepo.cpp
 * dnelson
 */
 
#include "transform.h"
#include "spectrum.h"
#include "volume.h"
#include "transfer.h"

#include "arepo.h"

#ifdef ENABLE_AREPO

// check for required Arepo compilation options

#ifndef VORONOI
#error ERROR. Missing required Arepo compilation option VORONOI.
#endif
#ifndef VORONOI_DYNAMIC_UPDATE
#error ERROR. Missing required Arepo compilation option VORONOI_DYNAMIC_UPDATE.
#endif

void Arepo::Init(int *argc, char*** argv)
{
		MPI_Init(argc, argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
		MPI_Comm_size(MPI_COMM_WORLD, &NTask);
		
		cout << "AREPO ENABLED. (NTask = " << NTask << " ThisTask = " << ThisTask << ")" << endl;
}

void Arepo::Cleanup()
{
		MPI_Finalize();
		close_logfiles();
}

bool Arepo::LoadSnapshot()
{
    IF_DEBUG(cout << "Arepo::LoadSnapshot(" << snapFilename << ")." << endl);
		
#ifndef DEBUG
		//freopen("/dev/null","w",stdout); //hide arepo stdout
#endif		
		
		// set startup options
		WriteMiscFiles = 0;
		RestartSnapNum = -1;
		RestartFlag    = SUNRISE_CODE;

		strcpy(ParameterFile,paramFilename.c_str());

		// call arepo: run setup
		begrun1();
		open_logfiles();
		
		// check snapshot exists
		if (ifstream(snapFilename.c_str())) {
				freopen("/dev/tty","w",stdout);
				cout << "Arepo::LoadSnapshot() ERROR! Exact snapshot [" << snapFilename << "] found (don't include extension)." << endl;
				terminate("1121");
		}
		
		string f1 = snapFilename + ".hdf5";
		string f2 = snapFilename + ".0.hdf5";
		if (!ifstream(f1.c_str()) && !ifstream(f2.c_str())) {
				freopen("/dev/tty","w",stdout);
				cout << "Arepo::LoadSnapshot() ERROR! Neither [" << f1 << "] nor [" << f2 << "] found!" << endl;
				terminate("1140");
		}
		
		// load snapshot (GAS ONLY)
		read_ic(snapFilename.c_str(), 0x01);
		
		// call arepo: read snapshot, allocate storage for tree, 
		//             initialize particle data, domain decomposition, initial HSML
  	if (init() != SUNRISE_CODE) {
				cout << "Arepo::LoadSnapshot() ERROR: Arepo did not return successfully." << endl;
				return false;
		}

#ifndef DEBUG
		//TODO: switch between these automatically
		//string fn = Config.imageFile + string(".out.txt");
		//freopen("/dev/tty","w",stdout); //return stdout to terminal (test only)
		//freopen(fn.c_str(),"a",stdout); //return stdout to a file (LSF)
#endif

		if (Config.verbose) {
				cout << endl << "Arepo Init Finished, Memory Report:" << endl;
				dump_memory_table();
		}
		
		return true;
}


ArepoMesh::ArepoMesh(const TransferFunction *tf)
{
		IF_DEBUG(cout << "ArepoMesh() constructor." << endl);
		
		// transfer function and sampling setup
		transferFunction = tf;
		viStepSize       = Config.viStepSize;
		
		//sampleWt = 0.2f; //All.BoxSize / pow(NumGas,0.333);
		sampleWt = 0.001f;
		
		if (viStepSize)
			sampleWt *= viStepSize;
		
		// set pointers into Arepo data structures
		T   = &Mesh;
		DP  = T->DP;
		DT  = T->DT;
		DTC = T->DTC;
		DTF = T->DTF;
		VF  = T->VF;
		//DC
		
		Ndp = T->Ndp;
		Ndt = T->Ndt;
		Nvf = T->Nvf;
		
		if (Config.verbose)
				cout << "[" << ThisTask << "] ArepoMesh: Ndp = " << Ndp << " Ndt = " << Ndt << " Nvf = " << Nvf 
						 << " NumGas = " << NumGas << " NumPart = " << NumPart << endl << endl;
		
		// boxsize
		extent = BBox(Point(0.0,0.0,0.0),Point(All.BoxSize,All.BoxSize,All.BoxSize));
		
		IF_DEBUG(extent.print(" ArepoMesh extent "));
		
		// preprocessing
		ArepoMesh::ComputeQuantityBounds();
		ArepoMesh::CalculateMidpoints();
		//ArepoMesh::LimitCellDensities();
		
		ArepoMesh::setupAuxMeshes();
		ArepoMesh::precomputeTetraGrads();
		
#ifdef NNI_WATSON_SAMBRIDGE
		// free volume accumulation array
		DP_vols = new float[Ndp];
#endif
		
		// TODO: temp units
		unitConversions[TF_VAL_DENS]   = All.UnitDensity_in_cgs / MSUN_PER_PC3_IN_CGS;
		unitConversions[TF_VAL_UTHERM] = All.UnitEnergy_in_cgs;
		
		IF_DEBUG(cout << "unitConv[dens]   = " << unitConversions[TF_VAL_DENS] << endl);
		IF_DEBUG(cout << "unitConv[utherm] = " << unitConversions[TF_VAL_UTHERM] << endl);
}

ArepoMesh::~ArepoMesh()
{				
#ifdef NATURAL_NEIGHBOR_INTERP
		// free aux meshes
		int numMeshes = Config.nTasks;
		
		for(int k = numMeshes-1; k == 0; k--)
		{
			myfree(AuxMeshes[k].DTF);
			myfree(AuxMeshes[k].DTC);
			AuxMeshes[k].DTC = NULL;
			myfree(AuxMeshes[k].DT);
			myfree(AuxMeshes[k].DP - 5);
			myfree(AuxMeshes[k].VF);
		}
		
		delete AuxMeshes;
#endif
#ifdef DTFE_INTERP
		// free delaunay tetra gradients
		delete DT_grad;
#endif
#ifdef NNI_WATSON_SAMBRIDGE
		// free volume accumulation array
		delete DP_vols;
#endif
}

void ArepoMesh::LocateEntryCellBrute(const Ray &ray)
{
		// note: using the brute force search is O(N_rays * NumPart) - not good
		Point hitbox  = ray(ray.min_t);
		
		int DP_ID = -1;
		double minDist = MAX_REAL_NUMBER;
		
		for (int i=0; i < Ndp; i++) {
				double dist = sqrt( (hitbox.x-DP[i].x)*(hitbox.x-DP[i].x) + 
													  (hitbox.y-DP[i].y)*(hitbox.y-DP[i].y) + 
														(hitbox.z-DP[i].z)*(hitbox.z-DP[i].z) );
				if (dist < minDist) {
						minDist = dist;
						DP_ID = i;
				}
		}

#ifdef DEBUG
		if (DP_ID >= NumGas)
				cout << " LocateEntryCellBrute is [Local Ghost] DP_ID = " << DP_ID << " (NumGas = " << NumGas << ")" << endl;
#endif
		
		IF_DEBUG(cout << " brute DP_ID = " << DP_ID << " DP[DP_ID].index = " << DP[DP_ID].index 
									<< " dist = " << minDist << " (DP.x = " << DP[DP_ID].x << " DP.y = " 
									<< DP[DP_ID].y << " DP.z = " << DP[DP_ID].z << ")" << endl);
		
		ray.index = DP_ID;
		ray.task  = DP[DP_ID].task;
		
		// verify task assignment
		if (ray.task < 0 || ray.task >= NTask) {
				cout << "ERROR! ray has bad task=" << ray.task << endl;
				terminate("1115");
		}
		
}

// set ray.min_t and ray.max_t bounds first

void ArepoMesh::LocateEntryCell(const Ray &ray, int *prevEntryCell)
{
		Point hitbox  = ray(ray.min_t);
		
#ifdef DEBUG
		Point exitbox = ray(ray.max_t);
		
		cout << " ray starts at x = " << hitbox.x << " y = " << hitbox.y << " z = " << hitbox.z << endl;
		cout << " ray ends at   x = " << exitbox.x << " y = " << exitbox.y << " z = " << exitbox.z << endl;
#endif		
		
		// use peanokey to find domain and task for ray
		if (ray.task == -1) {
				// TODO
				ray.task = 0;
		}
		// TODO: exchange
		
		// use tree to find nearest gas particle (local only)
		double mindist, mindist2;
		
#ifndef USE_AREPO_TREEFIND_FUNC
		const int dp_min = ArepoMesh::FindNearestGasParticle(hitbox, *prevEntryCell, &mindist);
		*prevEntryCell = dp_min;
#endif
		
#ifdef USE_AREPO_TREEFIND_FUNC
		terminate("6233");
		mindist = -1;
		double pos[3];
		pos[0] = hitbox.x;
		pos[1] = hitbox.y;
		pos[2] = hitbox.z;
		const int dp_min = 0; //= ngb_treefind_nearest_local(pos,-1,&mindist,0,&found_global);
#endif
	 
		IF_DEBUG(cout << " dp_min = " << dp_min
									<< " dist = " << mindist << " (x = " << P[dp_min].Pos[0] << " y = "
									<< P[dp_min].Pos[1] << " z = " << P[dp_min].Pos[2] << ")" << endl); 
	 
		// refine nearest point search to account for local ghosts
		int count   = 0;      // iterations marching through mesh
		int dp_old  = dp_min; // where we are
		int dp_oldi = dp_min; // candidates for closer points
		int dp_2ago = dp_min; // where we are coming from
		int dp_new;           // best candidate, where we are headed
		
		while (true)
		{
				// if any neighbors are closer, use them instead
				Vector celldist(hitbox.x - DP[dp_oldi].x,
											  hitbox.y - DP[dp_oldi].y,
												hitbox.z - DP[dp_oldi].z);
				
				mindist2 = celldist.LengthSquared();
				
				const int start_edge = midpoint_idx[dp_oldi].first;
				const int num_edges  = midpoint_idx[dp_oldi].second;
				
				IF_DEBUG(cout << " checking start_edge = " << start_edge << " num_edges = " << num_edges << endl);
				
				// search over all edges of this point
				for (int i=0; i < num_edges; i++)
				{
						const int dp_neighbor = opposite_points[start_edge + i];
						
						IF_DEBUG(cout << " iter i=" << i << " dp_neighbor = " << dp_neighbor << endl);
						
						// find distance to neighbor across this edge
						Point pos_neighbor(DP[dp_neighbor].x,DP[dp_neighbor].y,DP[dp_neighbor].z);
						
						celldist = hitbox - pos_neighbor;
						
						const float dist2 = celldist.LengthSquared();
						IF_DEBUG(cout << " dist2=" << dist2 << " mindist2=" << mindist2 << endl);
						
						if (dist2 < mindist2) {
								IF_DEBUG(cout << "  new closest DP_id = " << dp_neighbor << endl);
								mindist2 = dist2;
								dp_oldi = dp_neighbor;
						}
				}
				
				dp_new  = dp_oldi; // set best candidate, where we are going
		
				// prevent infinite loop, where we are headed where we just were
				if (count > 0 && dp_new == dp_2ago) {
						cout << "WARNING: LocateEntryCell refine bounce " << dp_2ago << " " << dp_new << endl;
						cout << " tree found dp_min=" << dp_min << " x = " << DP[dp_min].x << " y = " << 
										DP[dp_min].y << " z = " << DP[dp_min].z << endl;
						cout << " refine ended on dp_new=" << dp_new << " x = " << DP[dp_new].x << " y = " << 
										DP[dp_new].y << " z = " << DP[dp_new].z << endl;
						cout << " ray (hunting for, hitbox) at x = " << hitbox.x << " y = " << hitbox.y <<
										" z = " << hitbox.z << endl;
						continue;
				}
		
				dp_2ago = dp_old; // set where we are leaving
		
				// in closest if we didn't find any closer
				if (dp_new == dp_old) {
						IF_DEBUG(cout << " dp_new == dp_old = " << dp_new << " (in closest, entry search done)" << endl);
						break;
				}
				
				// not yet in closest, repeat search over edges for next closest cell
				IF_DEBUG(cout << " not yet in closest, moving to dp_new=" << dp_new << endl);
				dp_old = dp_new;
				count++;
				
				if ( count > 100 ) {
						cout << "Error: Refine treesearch hit iter=100." << endl;
						cout << " tree found dp_min=" << dp_min << " x = " << DP[dp_min].x << " y = " << 
						        DP[dp_min].y << " z = " << DP[dp_min].z << endl;
					  cout << " refine ended on dp_new=" << dp_new << " x = " << DP[dp_new].x << " y = " << 
						        DP[dp_new].y << " z = " << DP[dp_new].z << endl;
						cout << " ray (hunting for, hitbox) at x = " << hitbox.x << " y = " << hitbox.y <<
						        " z = " << hitbox.z << endl;
						terminate("1139");
				}
		}
		
		if (count > 20) {
				cout << "WARNING: LocateEntryCell iterated [" << count << "] times." << endl;
				cout << " tree found dp_min=" << dp_min << " x = " << DP[dp_min].x << " y = " << 
								DP[dp_min].y << " z = " << DP[dp_min].z << endl;
				cout << " refine ended on dp_new=" << dp_new << " x = " << DP[dp_new].x << " y = " << 
								DP[dp_new].y << " z = " << DP[dp_new].z << endl;
				cout << " ray (hunting for, hitbox) at x = " << hitbox.x << " y = " << hitbox.y <<
								" z = " << hitbox.z << endl;
		}
		
		// if we did not finish in a primary cell, check that we iterated over at least one neighbor
		if (dp_new >= NumGas && !count) {
				cout << "ERROR: Refined entry tree search ended in ghost but count=0" << endl;
				terminate("1107");
		}
		
		ray.index = dp_new;
		ray.task  = DP[dp_new].task;
		
		// verify task assignment
		if (ray.task < 0 || ray.task >= NTask) {
				cout << "ERROR! ray has bad task=" << ray.task << endl;
				terminate("1115");
		}
		
}

void ArepoMesh::VerifyPointInCell(int dp, Point &pos)
{
		Vector celldist(pos.x - DP[dp].x, pos.y - DP[dp].y, pos.z - DP[dp].z);
		
		double dist2point = celldist.LengthSquared();
		int dpid_min = dp;
		
#ifdef OLD_VERIFY
		// check neighbors for a closer tetra midpoint
		const int start_edge = midpoint_idx[dp].first;
		const int n_edges    = midpoint_idx[dp].second;
		
		for (int i=0; i < n_edges; i++)
		{
				const int dp_neighbor = opposite_points[start_edge + i];
				
				celldist = Vector(pos.x - DP[dp_neighbor].x,
				                  pos.y - DP[dp_neighbor].y,
													pos.z - DP[dp_neighbor].z);
				
				if (celldist.LengthSquared() < dist2point - INSIDE_EPS)
					dpid_min = dp_neighbor;
		}
#else
		// check all other DP for a closer tetra midpoint

		for (int i=0; i < Ndt; i++)
		{
			if ( DP[i].index >= NumGas )
			  continue; // skip ghosts
		
			celldist = Vector(pos.x - DP[i].x, pos.y - DP[i].y, pos.z - DP[i].z);
			
			if(celldist.LengthSquared() < dist2point - INSIDE_EPS)
				dpid_min = i;
		}
#endif
	
		// check
		if (dpid_min != dp) {
				cout << "VerifyPointInCell FAILED! dp = " << dp << " pos.x = " << pos.x << " pos.y = " << pos.y << " pos.z = " << pos.z << " (dist2point = " << dist2point << " DP.x = " << DP[dp].x << " DP.y = " << DP[dp].y << " DP.z = " << DP[dp].z << ")" << endl;
				
				Vector celldist(pos.x - DP[dpid_min].x, pos.y - DP[dpid_min].y, pos.z - DP[dpid_min].z);
				double dist2point = celldist.LengthSquared();
				
				cout << " dpid_min = " << dpid_min << " DP.x = " << DP[dpid_min].x << " DP.y = " << DP[dpid_min].y << " DP.z = " << DP[dpid_min].z << "(dist2point = " << dist2point << ")" << endl;
				terminate("1129");
		}
		
		IF_DEBUG(cout << "VerifyPointInCell Passed! dp = " << dp << " pos.x = " << pos.x << " pos.y = " << pos.y << " pos.z = " << pos.z << endl);

}

void ArepoMesh::LocateEntryTetra(const Ray &ray, int *prevEntryTetra)
{
		point p;
		int tt = 0, ret, moves;
		
		Point hitbox = ray(ray.min_t);
		
		// use get_tetra function with guess
		p.x = hitbox.x;
		p.y = hitbox.y;
		p.z = hitbox.z;
		set_integers_for_pointer(&p);
		
		tt = get_tetra(T, &p, &moves, *prevEntryTetra, &ret, &ret);
	 
#ifdef DEBUG
		int next_tetra = -1;
		cout << " TETRA tt = " << tt	<< " moves = " << moves << " ret = " << ret << endl;
		
		ret = -1;
		int check = InTetra(T, tt, &p, &ret, &next_tetra);
		cout << " check = " << check << " edgeface_nr = " << ret << " next_tetra = " << next_tetra << endl;
		if (check < 1) { // 2 or 3 ok (on face or edge)
			cout << "ERROR: Entry tetra failed." << endl;
			terminate("1164");
		}
#endif

		// sanity check
		if (tt < 0 || tt >= Ndt) {
				cout << "ERROR: Entry tetra search ended at tt = " << tt << endl;
				terminate("1107");
		}
		
		ray.tetra = tt;
		*prevEntryTetra = tt;
}

int ArepoMesh::FindNearestGasParticle(Point &pt, int guess, double *mindist)
{
	// based on ngbtree_walk.c:ngb_treefind_variable() (no MPI)
	int node, nearest, p;
	struct NgbNODE *current;
	double dx, dy, dz, cur_mindist, cur_mindist_sq, xtmp, ytmp, ztmp;
	float search_min[3], search_max[3], search_max_Lsub[3], search_min_Ladd[3];

#ifdef DEBUG
	int count_indpart=0,count_intnode=0,count_extnode=0;
#endif
	
	// starting node
	node = Ngb_MaxPart;
	
	if (guess >= 0) {
		nearest = guess;
	} else {
		// pick random gas particle for guess of the min distance
		//nearest = floor(get_random_number(SelRnd++) * NumGas);
		nearest = (int)floor(NumGas/2.0);
	}
	
	dx = NGB_PERIODIC_LONG_X(P[nearest].Pos[0] - pt.x);
	dy = NGB_PERIODIC_LONG_Y(P[nearest].Pos[1] - pt.y);
	dz = NGB_PERIODIC_LONG_Z(P[nearest].Pos[2] - pt.z);

	cur_mindist_sq = dx * dx + dy * dy + dz * dz;
  cur_mindist = sqrt(cur_mindist_sq);
	
	// search bounds
	search_min[0] = pt.x - cur_mindist;
	search_min[1] = pt.y - cur_mindist;
	search_min[2] = pt.z - cur_mindist;
	search_max[0] = pt.x + cur_mindist;
	search_max[1] = pt.y + cur_mindist;
	search_max[2] = pt.z + cur_mindist;

	search_max_Lsub[0] = search_max[0] - boxSize_X;
	search_max_Lsub[1] = search_max[1] - boxSize_Y;
	search_max_Lsub[2] = search_max[2] - boxSize_Z;

	search_min_Ladd[0] = search_min[0] + boxSize_X;
	search_min_Ladd[1] = search_min[1] + boxSize_Y;
	search_min_Ladd[2] = search_min[2] + boxSize_Z;
	

	while(node >= 0)
	{
		if(node < Ngb_MaxPart)  // single particle
		{
			IF_DEBUG(count_indpart++);
			p = node;
			node = Ngb_Nextnode[node];

			if(P[p].Type > 0) // not gas particle
				continue;
			
			dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - pt.x);
			if(dx > cur_mindist)
				continue;
				
			dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - pt.y);
			if(dy > cur_mindist)
				continue;
				
			dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - pt.z);
			if(dz > cur_mindist)
				continue;
	    
			double curdist2 = dx * dx + dy * dy + dz * dz;
			if(curdist2 > cur_mindist_sq)
				continue;
				
			cur_mindist_sq = curdist2;
			nearest = p;
		}
		else if(node < Ngb_MaxPart + Ngb_MaxNodes) // internal node
		{
			IF_DEBUG(count_intnode++);
			current = &Ngb_Nodes[node];

			// in case the node can be discarded
			node = current->u.d.sibling;

			// next check against bound
			if(search_min[0] > current->u.d.range_max[0] && search_max_Lsub[0] < current->u.d.range_min[0])
				continue;
			if(search_min_Ladd[0] > current->u.d.range_max[0] && search_max[0] < current->u.d.range_min[0])
				continue;

			if(search_min[1] > current->u.d.range_max[1] && search_max_Lsub[1] < current->u.d.range_min[1])
				continue;
			if(search_min_Ladd[1] > current->u.d.range_max[1] && search_max[1] < current->u.d.range_min[1])
				continue;

			if(search_min[2] > current->u.d.range_max[2] && search_max_Lsub[2] < current->u.d.range_min[2])
				continue;
			if(search_min_Ladd[2] > current->u.d.range_max[2] && search_max[2] < current->u.d.range_min[2])
				continue;

			node = current->u.d.nextnode; // need to open the node
		}
		else
		{
			IF_DEBUG(count_extnode++);
			node = Ngb_Nextnode[node - Ngb_MaxNodes];
			continue;
		}
	}

	if (nearest < 0 || nearest > NumGas) {
		cout << "ERROR: FindNearestGasParticle nearest=" << nearest << " out of bounds." << endl;
		terminate("1118");
	}
	
#ifdef DEBUG
	cout << "FindNearestGasParticle(): start nearest=" << nearest << " x = " << P[nearest].Pos[0] << " y = " << 
			P[nearest].Pos[1] << " z = " << P[nearest].Pos[2] << " cur_mindist = " << cur_mindist << endl;
	cout << " search_min = " << search_min[0] << " " << search_min[1] << " " << search_min[2] << endl;
	cout << " search_max = " << search_max[0] << " " << search_max[1] << " " << search_max[2] << endl;
	cout << " search_max_Lsub = " << search_max_Lsub[0] << " " << search_max_Lsub[1] << " " << search_max_Lsub[2] << endl;
	cout << " search_min_Ladd = " << search_min_Ladd[0] << " " << search_min_Ladd[1] << " " << search_min_Ladd[2] << endl;
	cout << " count_indpart = " << count_indpart << " count_intnode = " << count_intnode << " count_extnode = " << count_extnode << endl;
#endif	
	
	*mindist = sqrt(cur_mindist_sq);
	
	return nearest;
}

// get primary hydro ID - handle local ghosts
inline int ArepoMesh::getSphPID(int dp_id)
{
		int SphP_ID = -1;
		
		if (dp_id >= 0 && dp_id < NumGas)
				SphP_ID = dp_id;
		else if (dp_id >= NumGas)
				SphP_ID = dp_id - NumGas;
		
		if (SphP_ID < 0)
				terminate("1135");
				
		return SphP_ID;
}

void ArepoMesh::checkCurCellTF(bool *addFlag, int SphP_ID, float *vals)
{	
	// check if TF evaluates to zero at this cell midpoint
/*
	if (!transferFunction->InRange(vals)) {
			// if next cell in our path also evaluates TF to zero, skip this current cell, otherwise sample this
			// current cell to make sure we don't miss the ramp up of our TF
			// TODO: this is maybe still a good idea, but we need to check all neighboring
			//       cells, not just the one in the ray direction
			if (qmin != -1) {
					float vals_next[TF_NUM_VALS];
					vals_next[TF_VAL_DENS]        = (float) SphP[qmin].Density;
					vals_next[TF_VAL_UTHERM]      = (float) SphP[qmin].Utherm;
					
					if (!transferFunction->InRange(vals_next))
							addFlag = false;
			}
	}
*/		
}

void ArepoMesh::locateCurrentTetra(const Ray &ray, Vector &pt)
{
	// check degenerate point in R3, immediate return, otherwise we will terminate get_tetra
	// with "strange zero count" since we are on a vertex (3 faces simultaneously)
	if (fabs(pt.x - DP[ray.index].x) <= INSIDE_EPS && 
	    fabs(pt.y - DP[ray.index].y) <= INSIDE_EPS &&
			fabs(pt.z - DP[ray.index].z) <= INSIDE_EPS)
		return;
			
	// in making this substep, did we cross a delaunay tetra?
	point pp;
	pp.x = pt.x;
	pp.y = pt.y;
	pp.z = pt.z;
	set_integers_for_pointer(&pp);
			
	int next_tetra = 0,flag,edgeface_nr,moves;
			
	// check intersections of line defined by (pp,pend) with the faces of DT[ray.tetra]
	next_tetra = get_tetra(T, &pp, &moves, ray.tetra, &flag, &edgeface_nr);
	
//#ifdef DEBUG
	if ( next_tetra != ray.tetra )
	{
		//cout << "  TETRA ADVANCE old = " << ray.tetra << " new = " << next_tetra << endl;
				
		int ret,next_tetra2 = -1;
		int test = InTetra(T, next_tetra, &pp, &ret, &next_tetra2);
		if(test == 0)
		{
			cout << "  TETRA: ERROR, NOT INSIDE [" << test << "], wanted: " << next_tetra2 << endl;
			
			// just brute force figure out which tetra it really is in then
			int j = 0;
			for ( j = 0; j < Ndt; j++ )
			{
				if ( DT[j].p[0] == -5 || DT[j].p[1] == -5 || DT[j].p[2] == -5 || DT[j].p[3] == -5)
					continue;
							
				test = InTetra(T, j, &pp, &ret, &next_tetra2);
						
				if(test >= 1)
					cout << "    IN: [" << j << "] code: " << test << endl;
			}
							
			exit(20598);
		}
		
		//if(test > 1)
		//	cout << "  TETRA: WARNING: on face or edge [" << test << "]" << endl;
	}
//#endif

	ray.tetra = next_tetra;
}

bool ArepoMesh::AdvanceRayOneCellNew(const Ray &ray, float *t0, float *t1, 
																		 Spectrum &Lv, Spectrum &Tr, int taskNum)
{
		double min_t = MAX_REAL_NUMBER;
		int SphP_ID = getSphPID(DP[ray.index].index); // current primary cell index
		int qmin = -1; // next primary cell index
		
		// verify task
		if (ray.task != ThisTask) {
				cout << "[" << ThisTask << "] ERROR! ray.task = " << ray.task << " differs from ThisTask!" << endl;
				terminate("1138");
		}
	
	  Point pos = ray(ray.min_t);
	
#ifdef DEBUG_VERIFY_INCELL_EACH_STEP
		// verify ray is where we expect it
		ArepoMesh::VerifyPointInCell(ray.index,pos);
#endif

#ifdef DEBUG
		ray(*t0).print(" hitbox ");
		ray(*t1).print(" exitbox ");
#endif

		// DC connectivity

		double length; // find_next_voronoi_cell() return
		float p0[3];
		p0[0] = pos.x;
		p0[1] = pos.y;
		p0[2] = pos.z;
		
		double dir[3];
		dir[0] = ray.d[0];
		dir[1] = ray.d[1];
		dir[2] = ray.d[2];
		
		int prev_SphP_ID = getSphPID(DP[ray.prev_index].index);
		
		qmin = find_next_voronoi_cell(T, SphP_ID, p0, dir, prev_SphP_ID, &length); //&(pos[0])
		qmin = DC[qmin].index;
		min_t = ray.min_t + length;
		
		IF_DEBUG(cout << "  NEW intersection t = " << min_t << " setting new min_t, qmin (next DP) = " << qmin << endl);
		
		// alternative Sunrise connectivity
		int qmin_old = -1;
		int min_t_old = -1;
		Point hitbox  = ray(*t0);
		const Vector cellp(DP[ray.index].x,DP[ray.index].y,DP[ray.index].z);
		
		const pair<int,int> edge = midpoint_idx[ray.index];
		
		for (int i=edge.second-1; i >= 0; i--)
		{
			// skip face we arrived through, if any
			if (opposite_points[edge.first + i] == ray.prev_index && ray.prev_index != -1)
				continue;
			
			//IF_DEBUG(cout << " OLD checking face[" << i << "] midpoints.x = " << midpoints[edge.first + i].x
			//							<< " midpoints.y = " << midpoints[edge.first + i].y << " midpoints.z = " 
			//							<< midpoints[edge.first + i].z << " opposite_point = " 
			//							<< opposite_points[edge.first + i] << endl);
	
			// midpoint (c)
			const Vector midp(midpoints[edge.first + i].x - hitbox.x,
												midpoints[edge.first + i].y - hitbox.y,
												midpoints[edge.first + i].z - hitbox.z); //TODO: verify hitbox not hitcell
			
			// vector pointing to the outside, normal to a voronoi face of the cell (q)
			const Vector norm = midpoints[edge.first + i] - cellp;
	
			//IF_DEBUG(midp.print("  midp (midpoints-hitbox) "));
			//IF_DEBUG(norm.print("  norm (midpoints-cellp) "));
	
			// find intersection of ray with this face
			double dotprod1 = Dot( ray.d, norm ); //ddotq
			double dotprod2 = Dot( midp,  norm ); //cdotq
			
			// check if ray is aligned on face (e.g. backgroundgrid)
			if (dotprod1 == 0 && dotprod2 == 0)
				continue;
			
			if (dotprod1 > 0)
			{
				double t = dotprod2 / dotprod1; // infinite line/plane intersection test
				
				//IF_DEBUG(cout << "  OLD i[" << i << "] dotprod>0 and t = " << t << " (min=" << (ray.min_t-*t0) << ")" << endl);
				
				if (t > (ray.min_t-*t0) && t < min_t) {
					min_t_old = t;
					qmin_old = opposite_points[edge.first + i];
					IF_DEBUG(cout << "  OLD intersection t = " << min_t_old << " setting new min_t, qmin (next DP) = " << qmin_old << endl);
				}
			}
		}
		
		IF_DEBUG(cout << "  qmin = " << qmin << " qmin_old = " << qmin_old << " min_t = " << min_t << " min_t_old = " << min_t_old << endl);
		
#ifdef USE_SUNRISE_CONNECTIVITY
		qmin = qmin_old;
		min_t = min_t_old;
#endif
	
		// check if exiting box and failed to exit a face
		if (qmin == -1) {
			Point exitcell = ray(*t0 + min_t);
			
			if (!extent.Inside(exitcell)) { // && ray.index >= NumGas
				// set intersection with box face to allow for final contribution to ray
				IF_DEBUG(cout << " failed to intersect face, exitcell outside box, ok!" << endl);
				min_t = ray.max_t;
				
				// fake exit face
				qmin = 0;
			}
		}
		
		// disable adding any contribution from ghosts to rays
		bool addFlag = true;
		if (qmin != -1 && ray.index >= NumGas)
			addFlag = false;
		
		// pack cell-center values to test TF
		float vals[TF_NUM_VALS];
		
		vals[TF_VAL_DENS]        = (float) SphP[SphP_ID].Density;
		vals[TF_VAL_UTHERM]      = (float) SphP[SphP_ID].Utherm;
		vals[TF_VAL_PRES]        = (float) SphP[SphP_ID].Pressure;
		vals[TF_VAL_ENERGY]      = (float) SphP[SphP_ID].Energy;
		
		vals[TF_VAL_VEL_X]       = (float) P[SphP_ID].Vel[0];
		vals[TF_VAL_VEL_Y]       = (float) P[SphP_ID].Vel[1];
		vals[TF_VAL_VEL_Z]       = (float) P[SphP_ID].Vel[2];
		
		// if TF evaluates to zero at this cell midpoint, and at all neighboring cell midpoints, 
		// then we are guaranteed that nowhere in this cell will be important, so jump onwards
		checkCurCellTF(&addFlag, SphP_ID, vals);

		// check for proper exit point
		if (qmin != -1)
		{
			// clamp min_t to avoid integrating outside the box
			IF_DEBUG(cout << " min_t = " << min_t << " (t1=" << *t1 << " t0=" << *t0 << ")" << endl);
			min_t = Clamp(min_t,0.0,(*t1-*t0));
			
			if (addFlag)
			{
				// entry and exit points for this cell
				Point hitcell  = ray(ray.min_t);
				Point exitcell = ray(*t0 + min_t);

				IF_DEBUG(hitcell.print(" hcell "));
				IF_DEBUG(exitcell.print(" ecell "));
					
				// cell gradient information
				Vector sphCen(     SphP[SphP_ID].Center[0],   SphP[SphP_ID].Center[1],   SphP[SphP_ID].Center[2]);
				Vector sphDensGrad(SphP[SphP_ID].Grad.drho[0],SphP[SphP_ID].Grad.drho[1],SphP[SphP_ID].Grad.drho[2]);
															 
				const Vector norm = exitcell - hitcell;

				IF_DEBUG(norm.print(" norm "));
				IF_DEBUG(sphCen.print(" sphCen "));
					
				// compute total path length through cell
				double len = norm.Length();
					
				// find interpolated density at midpoint of line segment through voronoi cell
				Vector midpt(hitcell[0] + 0.5 * norm[0],hitcell[1] + 0.5 * norm[1],hitcell[2] + 0.5 * norm[2]);
				midpt -= sphCen;

				// optical depth: always use gradient for tau calculation (though apply as one value for whole cell)
				Spectrum stepTau(0.0);
				stepTau += transferFunction->sigma_t() * (SphP[SphP_ID].Density + Dot(sphDensGrad,midpt)) * len;
				//Tr *= Exp(-stepTau); // reduce transmittance for optical depth
					
				// TODO: sub-step length should be adaptive based on gradients

				// if not sub-stepping then set default
				if (!viStepSize)
					viStepSize = len;
						
				// sub-stepping: variable 
				//int nSamples = (int)ceilf(len_sample / viStepSize);
				//double fracstep = 1.0 / nSamples;
				//double halfstep = 0.5 / nSamples;
					
				// setup sub-stepping: strict in world space
				Vector prev_sample_pt(ray(ray.depth * viStepSize));
				Vector norm_sample(exitcell - prev_sample_pt);

				int nSamples = (int)floorf(norm_sample.Length() / viStepSize);
				norm_sample = Normalize(norm_sample);
					
				IF_DEBUG(prev_sample_pt.print(" prev_sample_pt "));
					
				IF_DEBUG(cout << " sub-stepping len = " << len << " nSamples = " << nSamples 
											<< " (step = " << len/nSamples << ")" << endl);
												
				for (int i = 0; i < nSamples; ++i)
				{
					// where are we inside the cell?
					Vector midpt(prev_sample_pt[0] + ((i+1)*viStepSize) * norm_sample[0],
											 prev_sample_pt[1] + ((i+1)*viStepSize) * norm_sample[1],
											 prev_sample_pt[2] + ((i+1)*viStepSize) * norm_sample[2]);
							
					IF_DEBUG(midpt.print("  substep midpt "));
					
#if defined(DTFE_INTERP) || defined(NNI_WATSON_SAMBRIDGE) || defined(NNI_LIANG_HALE)
					locateCurrentTetra(ray, midpt);
#endif
																							
					// subsample (replace fields in vals by interpolated values)
					int status = subSampleCell(SphP_ID, ray, midpt, &vals[0], taskNum);
					if(status) status = 0; // kill warning
							
#ifdef DEBUG
					double fracstep = 1.0 / nSamples;
					cout << "  segment[" << i << "] fractrange [" << (i*fracstep) << "," 
							<< (i*fracstep)+fracstep << "] rho = " << SphP[SphP_ID].Density
							<< " rho subSample = " << vals[TF_VAL_DENS] << endl;
#endif
							
					// apply TF to integrated (total) quantities (only appropriate for constant?)
					if (Config.projColDens) {
						terminate("1299"); // best check this
						vals[TF_VAL_DENS] *= len;
					}
					
					// compute emission-only source term using transfer function
					Lv += Tr * transferFunction->Lve(vals) * sampleWt;
							
					// update previous sample point marker
					ray.depth++;
				} // nSamples
			} //addFlag
			
			// update ray: transfer to next voronoi cell (possibly on different task)
			ray.task  = DP[qmin].task;
			ray.prev_index = ray.index;
			ray.index = qmin;
			ray.min_t = Clamp(min_t + *t0,ray.min_t,ray.max_t);
			
			IF_DEBUG(cout << " updated ray new task = " << ray.task << " index = " << ray.index 
										<< " min_t = " << ray.min_t << endl);
			
			if (fabs(ray.min_t - ray.max_t) <= INSIDE_EPS) {
					// apparently this ray is done?
					IF_DEBUG(cout << " min_t == t1 = " << *t1 << ", ray done." << endl);
					return false;
			}
		
		} else {
			// failed to intersect a face (only should happen if exiting box - no connectivity with big tetra)
				
			if (ray.min_t < ray.max_t - INSIDE_EPS) {
					// in primary cell or exitpoint inside box, either way this should not happen
					cout << "ERROR! Ray did not finish. min_t = " << ray.min_t << " max_t = " << ray.max_t << endl;
					cout << " DP[ray.index] pos: " << DP[ray.index].x << " " << DP[ray.index].y 
							 << " " << DP[ray.index].z << endl << endl;
					terminate("1130");
			}
		} // qmin
		
		return true;		
}

// for now just zero hydro quantities of primary cells that extend beyond the box
void ArepoMesh::LimitCellDensities()
{
		// loop over all tetras
		for (int i=0; i < Ndt; i++)
		{
				// skip those with initial points outside the box or connecting to DPinfinity
				if (DT[i].t[0] < 0 || DT[i].p[0] == DPinfinity || DT[i].p[1] == DPinfinity
				                   || DT[i].p[2] == DPinfinity || DT[i].p[3] == DPinfinity)
						continue;
						
				// circumsphere center
				Point dtc(DTC[i].cx,DTC[i].cy,DTC[i].cz);
				
				// loop over the 4 vertices
				for (int j=0; j < 4; j++)
				{
						// find the cell opposite this vertex
						const int dp = DT[i].p[j];
				
						int SphP_ID = getSphPID(DP[dp].index);
								
						// valid cell?
						if (DP[dp].index < NumGas && SphP_ID >= 0) {
								if (!extent.Inside(dtc)) {
										IF_DEBUG(cout << " Zeroing Density and Grad SphP_ID=" << SphP_ID << " dtc.x = " << dtc.x
										              << " dtc.y = " << dtc.y << " dtc.z = " << dtc.z << endl);
																	
										SphP[SphP_ID].Density = 0;
										SphP[SphP_ID].Grad.drho[0] = 0;
										SphP[SphP_ID].Grad.drho[1] = 0;
										SphP[SphP_ID].Grad.drho[2] = 0;
										//TODO: zero any other quantitfy used in a TF (e.g. utherm)
								}
						} // valid?
				} // vertices
		} //tetras

}

// construct the sunrise alternative connectivity

/*  Sets up the midpoints array which stores the midpoints of the face
    planes in a compact way so we don't have to chase a bunch of
    pointers all over the place to do the intersection tests with the
    face planes when finding the voronoi neighbors. The connections
    stored in our data structure differs from that in the Arepo DC
    array in that we map connections between distinct mesh points, not
    just primary cells. */
void ArepoMesh::CalculateMidpoints()
{
		IF_DEBUG(cout << "ArepoMesh::CalculateMidpoints()" << endl);

		// set up temporary multimap of SphP id -> DP id, to identify all local ghosts
		// associated with a particular SphP entry
		multimap<int,int> mm;
		typedef multimap<int,int>::iterator mmi;
		
		for (int i=0; i < Ndp; i++)
		{
				int SphP_ID = getSphPID(DP[i].index);
				
				mm.insert(make_pair(SphP_ID, i));
		}
		
		// set up mapping of DP id -> DP primary id
		for (int i=0; i < Ndp; i++)
		{
				int SphP_ID = getSphPID(DP[i].index);
		
				// cell has no hydro quantities -> map to -1
				if (SphP_ID < 0) {
						IF_DEBUG(cout << "WARNING: CM i=" << i << " SphP_ID (neg) = " << SphP_ID << endl);
						primary_cells.push_back(-1);
				}
				
				// loop over all DP indices that share this SphP cell
				pair<mmi,mmi> dp_indices(mm.equal_range(SphP_ID));
				
				if (dp_indices.second == dp_indices.first)
						terminate("1131");
						
				// search for primary cell and record to map
				while (dp_indices.first != dp_indices.second)
				{
						if (dp_indices.first->second >= 0 && dp_indices.first->second < NumGas) {
								primary_cells.push_back(dp_indices.first->second);
								break;
						}
						dp_indices.first++;
				}
		}
		
		// verify size
		if (primary_cells.size() != (unsigned int)Ndp)
				terminate("1132");
				
		// use VF array to generate connnections (reorganize it to index by point)
		multimap<int,int> conn;
		
		for (int i=0; i < Nvf; i++) {
				conn.insert(make_pair(VF[i].p1,VF[i].p2));
				conn.insert(make_pair(VF[i].p2,VF[i].p1));
		}
		
		// associate connections with cells
		for (int i=0; i < Ndp; i++)
		{
				//int SphP_ID = getSphPID(DP[i].index);
						
				const Vector cellp(DP[i].x,DP[i].y,DP[i].z);
				
				// find connections for this cell
				pair<mmi,mmi> dp_neighbors(conn.equal_range(i));
				
				if (dp_neighbors.first == dp_neighbors.second)
						terminate("1133");
						
				for (; dp_neighbors.first != dp_neighbors.second; dp_neighbors.first++)
				{
						const int dp_neighbor = dp_neighbors.first->second;
				
						int SphP_ID_n = -1; // don't use getSphPID()
						
						if (DP[dp_neighbor].index >= 0 && DP[dp_neighbor].index < NumGas)
								SphP_ID_n = DP[dp_neighbor].index;
						else if (DP[dp_neighbor].index >= NumGas) {
								SphP_ID_n = DP[dp_neighbor].index - NumGas;
						}
								
						// skip invalid neighbors
						if (SphP_ID_n < 0)
								continue;
								
						const Vector midp( 0.5 * (cellp.x + DP[dp_neighbor].x),
															 0.5 * (cellp.y + DP[dp_neighbor].y),
															 0.5 * (cellp.z + DP[dp_neighbor].z) );
																										 
						// add connection to the connectivity map
						midpoints.push_back(midp);
						opposite_points.push_back(dp_neighbor);
				}
				
				// all connections for this DP cell done, update the number and position in the midpoint index vector
				const int start_pos = 
					midpoint_idx.empty() ? 0 :
					(midpoint_idx.back().first + midpoint_idx.back().second);
				midpoint_idx.push_back(make_pair(start_pos, midpoints.size()-start_pos));
		}

}

void ArepoMesh::ComputeQuantityBounds()
{
		float pmax  = 0.0;
		float pmin  = INFINITY;
		float pmean = 0.0;
		
		float umax = 0.0;
		float umin = INFINITY;
		float umean = 0.0;
		
		for (int i=0; i < NumGas; i++) {
				if (SphP[i].Density > pmax)
						pmax = SphP[i].Density;
				if (SphP[i].Density < pmin)
						pmin = SphP[i].Density;
				pmean += SphP[i].Density;
				
				if (SphP[i].Utherm > umax)
						umax = SphP[i].Utherm;
				if (SphP[i].Utherm < umin)
						umin = SphP[i].Utherm;
				umean += SphP[i].Utherm;
		}
		
		valBounds[TF_VAL_DENS*3 + 0] = pmin;
		valBounds[TF_VAL_DENS*3 + 1] = pmax;
		valBounds[TF_VAL_DENS*3 + 2] = pmean / NumGas;
		
		valBounds[TF_VAL_UTHERM*3 + 0] = umin;
		valBounds[TF_VAL_UTHERM*3 + 1] = umax;
		valBounds[TF_VAL_UTHERM*3 + 2] = umean / NumGas;
		
		cout << " Density min = " << valBounds[TF_VAL_DENS*3 + 0] 
									<< " max = " << valBounds[TF_VAL_DENS*3 + 1] 
									<< " mean = " << valBounds[TF_VAL_DENS*3 + 2] << endl;
									
		cout << " Utherm  min = " << valBounds[TF_VAL_UTHERM*3 + 0] 
									<< " max = " << valBounds[TF_VAL_UTHERM*3 + 1] 
									<< " mean = " << valBounds[TF_VAL_UTHERM*3 + 2] << endl;
			
/*			
		for (int i = 0; i < NumGas; i++) {
				//SphP[i].Density /= valBounds[TF_VAL_DENS*3 + 0];
				//SphP[i].Utherm  /= valBounds[TF_VAL_UTHERM*3 + 0];
				//SphP[i].Density = log(SphP[i].Density);
				//SphP[i].Utherm  = log(SphP[i].Utherm);
				SphP[i].Density = 1.0;
				SphP[i].Utherm  = 1.0;
		}
	*/	

	/*
		float invMaxMinusMin = 1.0f / (valBounds[TF_VAL_DENS*3 + 1] - valBounds[TF_VAL_DENS*3 + 0]);
	
		for (int i = 0; i < NumGas; i++) {
				SphP[i].Density = (SphP[i].Density - valBounds[TF_VAL_DENS*3 + 0]) * invMaxMinusMin;
		}	
		
		*/
}

int ArepoMesh::ComputeVoronoiEdges()
{
		IF_DEBUG(cout << "ArepoMesh::ComputeVoronoiEdges()" << endl);

		// geometric conventions (voronoi.h)
		const int edge_start[6]     = { 0, 0, 0, 1, 1, 2 };
		const int edge_end[6]       = { 1, 2, 3, 2, 3, 3 };
		const int edge_opposite[6]  = { 3, 1, 2, 3, 0, 1 };
		const int edge_nexttetra[6] = { 2, 3, 1, 0, 2, 0 };		
		
		tetra *prev, *next;
		int i,j,k,l,m,ii,jj,kk,ll,tt,next_tt;
		int dp1,dp2,edge_nr,bit,nr_next,count;
		
		vertexList.reserve(2*Nvf);
		numVertices.reserve(Nvf);
		vertexOffset.reserve(Nvf);
		
		Edge_visited = new unsigned char[Ndt];
		
		// zero
		for(i = 0; i < Ndt; i++)
			Edge_visited[i] = 0;

		// loop over all local tetra
		for(tt = 0; tt < Ndt; tt++)
		{
		
			if (Mesh.DT[tt].t[0] < 0) // skip deleted tetras
				continue;
	
			bit  = 1;
			edge_nr = 0;
			
			// loop over all edges of this tetra
			while( Edge_visited[tt] != EDGE_ALL ) {
			
				if( (Edge_visited[tt] & bit) != 0 ) {
					bit <<= 1;
					edge_nr++;
					continue;
				}

				tetra *t = &DT[tt];

				// edge-point relation
				i = edge_start[edge_nr];
				j = edge_end[edge_nr];
				k = edge_opposite[edge_nr];
				l = edge_nexttetra[edge_nr];
			
				// mark edge as visited
				Edge_visited[tt] |= (1 << edge_nr);

				// delaunay points on both side of face
				dp1 = t->p[i];
				dp2 = t->p[j];

				// skip large tetra
				if(dp1 < 0 || dp2 < 0) {
					bit <<= 1;
					edge_nr++;
					continue;
				}

				// skip ghost points (both local and foreign)
				if((DP[dp1].task != ThisTask || DP[dp1].index < 0 || DP[dp1].index >= NumGas) &&
					 (DP[dp2].task != ThisTask || DP[dp2].index < 0 || DP[dp2].index >= NumGas)) {
					bit <<= 1;
					edge_nr++;
					continue;
				}

				// count number of face vertices
				count = 0;
				prev = t;

				do
				{
					count++;
					next_tt = prev->t[l];
					next = &DT[next_tt];

					for(m = 0, ll = ii = jj = -1; m < 4; m++) {
						if(next->p[m] == prev->p[k])
							ll = m;
						if(next->p[m] == prev->p[i])
							ii = m;
						if(next->p[m] == prev->p[j])
							jj = m;
					}

					if(ll < 0 || ii < 0 || jj < 0)
						terminate("inconsistency");

					kk = 6 - (ll + ii + jj);
					i = ii;
					l = ll;
					j = jj;
					k = kk;

					prev = next;
				}
				while(next != t);

				count++;

				// add count of vertices for this face to Nvertices and first vertex tetra index to VertexList				
				numVertices.push_back(count);
				vertexList.push_back(tt);
				
				IF_DEBUG(cout << " face i=" << numVertices.size() << " have [" << count << "] vertices" << endl);
				
				// add subsequent tetra indices for the other vertices of this voronoi face
				count = 0;
				prev = t;
				
				do
				{
					count++;
					next_tt = prev->t[l];
					next = &DT[next_tt];
		
					vertexList.push_back(next_tt);
						
					IF_DEBUG(cout << "  adding to face i=" << numVertices.size() << " VertexList[" << vertexList.size() << "] = " << vertexList.back() << endl);
		
					for(m = 0, ll = ii = jj = -1; m < 4; m++)
					{
						if(next->p[m] == prev->p[k])
							ll = m;
						if(next->p[m] == prev->p[i])
							ii = m;
						if(next->p[m] == prev->p[j])
							jj = m;
					}

					if(ll < 0 || ii < 0 || jj < 0)
						terminate("inconsistency");

					kk = 6 - (ll + ii + jj);

					// flag edge
					for(nr_next = 0; nr_next < 6; nr_next++) {
						if((edge_start[nr_next] == ii && edge_end[nr_next] == jj) || (edge_start[nr_next] == jj && edge_end[nr_next] == ii))
						{
							if((Edge_visited[next_tt] & (1 << nr_next)) && next != t)
								terminate("inconsistency");

							Edge_visited[next_tt] |= (1 << nr_next);
							break;
						}
					}

					i = ii;
					l = ll;
					j = jj;
					k = kk;

					prev = next;

				}
				while(next != t);

				bit <<= 1;
				edge_nr++;
		
			} // edges
		} // tt		
		
		// create offset table
		vertexOffset.push_back(0);
		
		for(size_t i = 1; i < numVertices.size(); i++) {
			vertexOffset.push_back(vertexOffset[i - 1] + numVertices[i - 1]);
			cout << "[" << i << "] numVert=" << numVertices[i-1] << " offset=" << vertexOffset.back() << endl;
		}
		
		delete Edge_visited;
		
		return numVertices.size();
		
}

#ifdef DUMP_VORONOI_MESH
void ArepoMesh::OutputMesh()
{
	char buf[500];
	sprintf(buf,"voronoi_mesh_0");
	write_voronoi_mesh(&Mesh,buf,0,0);
	cout << "MESH WRITTEN." << endl;
}
#endif

void ArepoMesh::DumpMesh()
{
		cout << endl << "Delaunay Points [" << Ndp << "]:" << endl;
		for (int i=0; i < Ndp; i++) {
				cout << setw(3) << i << " x = " << DP[i].x 
														 << " y = " << DP[i].y 
														 << " z = " << DP[i].z 
														 << " xx = " << DP[i].xx 
														 << " yy = " << DP[i].yy 
														 << " zz = " << DP[i].zz << endl
														 << "    ID = " << DP[i].ID 
														 << " task = " << DP[i].task 
														 << " index = " << DP[i].index 
														 << " ix = " << DP[i].ix 
														 << " iy = " << DP[i].iy 
														 << " iz = " << DP[i].iz << endl;
		}
		
		cout << endl << "SphP Hydro [" << NumGas << "]:" << endl;
		for (int i=0; i < NumGas; i++) {
				cout << setw(3) << i << " dens = " << SphP[i].Density
														 << " pres = " << SphP[i].Pressure
														 << " uthm = " << SphP[i].Utherm
														 << " energy = " << SphP[i].Energy
														 << " p[0] = " << SphP[i].Momentum[0]
														 << " p[1] = " << SphP[i].Momentum[1]
														 << " p[2] = " << SphP[i].Momentum[2]
														 << " vol = " << SphP[i].Volume
														 << " oldmass = " << SphP[i].OldMass << endl;		
		}

		cout << endl << "Delaunay Tetra [" << Ndt << "] [DIMS=" << DIMS << "]:" << endl;
		for (int i=0; i < Ndt; i++) {
				cout << setw(3) << i;
				for (int j=0; j < DIMS+1; j++)
				  cout << " p[" << j << "] = " << setw(2) << DT[i].p[j];
				for (int j=0; j < DIMS+1; j++)
				  cout << " t[" << j << "] = " << setw(2) << DT[i].t[j];
				for (int j=0; j < DIMS+1; j++)
				  cout << " s[" << j << "] = " << setw(1) << DT[i].s[j];
				cout << endl;
		}
		
		cout << endl << "Delaunay Circumcircle Centers:" << endl;
		for (int i=0; i < Ndt; i++) {
				cout << setw(3) << i << " cx = " << setw(8) << DTC[i].cx 
														 << " cy = " << setw(8) << DTC[i].cy 
														 << " cz = " << setw(8) << DTC[i].cz << endl;
		}
		
		//cout << endl << "Delaunay Tetra Faces:" << endl;
		//for (int i=0; i < Ndt; i++) {
		//		cout << setw(3) << i << " DTF = " << DTF[i] << endl;
		//}
		
		cout << endl << "Voronoi Faces [" << Nvf << "]:" << endl;
		for (int i=0; i < Nvf; i++) {
				cout << setw(3) << i << " p1 = " << setw(3) << VF[i].p1
				                     << " p2 = " << setw(3) << VF[i].p2
														 << " area = " << setw(12) << VF[i].area
														 << " cx = " << setw(10) << VF[i].cx
														 << " cy = " << setw(10) << VF[i].cy
														 << " cz = " << setw(10) << VF[i].cz << endl;
		
		}

    cout << endl << "Voronoi Connections (DC):" << endl;
		int dc,next;
		for (int i=0; i < NumGas; i++) {
				dc   = SphP[i].first_connection;
				cout << " SphP[" << setw(3) << i << "] DC.first = " << setw(2) << SphP[i].first_connection;
				
				do {
						next = DC[dc].next;
						cout << " " << " next = " << setw(2) << next;
						dc = next;
				} while (next != SphP[i].last_connection);
				
				cout << " DC.last = " << SphP[i].last_connection << endl;		
		}
		
		if (numVertices.size()) {
			cout << endl << "Voronoi Edges (NumGas=" << NumGas << "):" << endl;
			for (int i=0; i < NumGas; i++) {
					cout << setw(3) << i << " numVert = " << setw(2) << numVertices[i] << " vertexOffset = " 
							 << setw(3) << vertexOffset[i] << endl;
			}
		}
		
		cout << endl << "Primary_Cells and Midpoint_Idx (size=" << primary_cells.size() << "):" << endl;
		for (unsigned int i=0; i < primary_cells.size(); i++) {
				cout << "[" << setw(2) << i << "] primary id = " << primary_cells[i] 
				     << " edges start " << midpoint_idx[i].first << " num edges = " 
						 << midpoint_idx[i].second << endl;
		}
		
		cout << endl << "Midpoints and Opposite_Points (size=" << midpoints.size() << "):" << endl;
		for (unsigned int i=0; i < opposite_points.size(); i++) {
				cout << "[" << setw(2) << i << "] x = " << midpoints[i].x << " y = " << midpoints[i].y
				     << " z = " << midpoints[i].z << " opposite id = " << opposite_points[i] << endl;
		}
}

bool ArepoMesh::TetraEdges(const int i, vector<Line> *edges)
{
		Point pts[4];
		
		for (int j=0; j < DIMS+1; j++) {
				// if this tetra includes the "infinity" point do not draw
				//if (DT[i].p[j] == -5) {
				// if this tetra includes the inf point or any of the global bounding tetra, do not draw
				if (DT[i].p[j] < 0) {
						IF_DEBUG(cout << " edge[" << i << "] pt[" << j << "] is INFINITY, skipping." << endl);
						return false;
				}				
					
				pts[j] = Point(DP[DT[i].p[j]].x,
											 DP[DT[i].p[j]].y, 
											 DP[DT[i].p[j]].z);													
				IF_DEBUG(cout << " edge[" << i << "] pt[" << j << "] DP ind = " << DT[i].p[j] 
      			 << " x = " << pts[j].x << " y = " << pts[j].y << " z = " << pts[j].z << endl);
		}
			
		edges->push_back(Line(pts[0],pts[1])); // 1 - base one
		edges->push_back(Line(pts[1],pts[2])); // 2 - base two
		edges->push_back(Line(pts[2],pts[0])); // 3 - base three
		
		edges->push_back(Line(pts[0],pts[3])); // 4 - (0,top)
		edges->push_back(Line(pts[1],pts[3])); // 5 - (1,top)
		edges->push_back(Line(pts[2],pts[3])); // 6 - (2,top)
		
		return true;
}

bool ArepoMesh::VoronoiEdges(const int i_face, vector<Line> *edges)
{
		IF_DEBUG(cout << "VoronoiEdges(" << i_face << ") numVertices=" << numVertices[i_face]
									<< " vertexOffset = " << vertexOffset[i_face] << endl);
		
		if (numVertices[i_face] <= 0 || numVertices[i_face] < DIMS || i_face < 0 || i_face >= Nvf) {
				IF_DEBUG(cout << "WARNING: Nvert[" << i_face << "] empty, degenerate, or out of bounds." << endl);
				return false;
		}
		
		int s_ind = vertexList[vertexOffset[i_face]];
		int n_ind;
		
		Point prev = Point(DTC[s_ind].cx, DTC[s_ind].cy, DTC[s_ind].cz);
		Point next;

		// loop over remaining vertices (one extra at end with modulo to connect to first point)
		for (int i=1; i < numVertices[i_face]; i++) {	
				n_ind = vertexList[(vertexOffset[i_face] + i)];// % numVertices[i_face]];
				next  = Point(DTC[n_ind].cx, DTC[n_ind].cy, DTC[n_ind].cz);
				
				if (!extent.Inside(prev) || !extent.Inside(next)) {
						IF_DEBUG(cout << " VE[" << i << "] circumcircle center outside extent, skipping." << endl);
						continue;
				}
				
				edges->push_back(Line(prev,next));
				
				IF_DEBUG(cout << " edge[" << i_face << "," << i 
											<< "] prev.x = " << prev.x << " prev.y = " << prev.y << " prev.z = " << prev.z
											<<  " next.x = " << next.x << " next.y = " << next.y << " next.z = " << next.z << endl);
						 
				prev = next;
		}
		
		return true;
}

#endif //ENABLE_AREPO
