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
//#ifndef VORONOI_NEW_IMAGE
//#error ERROR. Missing required Arepo compilation option VORONOI_NEW_IMAGE.
//#endif
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
		RestartFlag    = SUNRISE;

		strcpy(ParameterFile,paramFilename.c_str());

		// call arepo: run setup
		begrun1();
    
		// check snapshot exists
		if (ifstream(snapFilename.c_str())) {
				freopen("/dev/tty","w",stdout);
				cout << "Arepo::LoadSnapshot() ERROR! Exact snapshot [" << snapFilename << "] found (don't include extension)." << endl;
				endrun(1121);
		}
		
		string f1 = snapFilename + ".hdf5";
		string f2 = snapFilename + ".0.hdf5";
		if (!ifstream(f1.c_str()) && !ifstream(f2.c_str())) {
				freopen("/dev/tty","w",stdout);
				cout << "Arepo::LoadSnapshot() ERROR! Neither [" << f1 << "] nor [" << f2 << "] found!" << endl;
				endrun(1140);
		}
		
		// load snapshot (GAS ONLY)
		read_ic(snapFilename.c_str(), 0x01);
		
		// call arepo: read snapshot, allocate storage for tree, 
		//             initialize particle data, domain decomposition, initial HSML
  	if (init() != SUNRISE) {
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
		
		// transfer function and transform
		transferFunction = tf;
		viStepSize       = Config.viStepSize;
		
		// set NULL
		Nedges       = NULL;
		EdgeList     = NULL;
		NedgesOffset = NULL;
		
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
						 << " N_gas = " << N_gas << " NumPart = " << NumPart << endl << endl;
		
		// boxsize
		extent = BBox(Point(0.0,0.0,0.0),Point(All.BoxSize,All.BoxSize,All.BoxSize));
		
		IF_DEBUG(extent.print(" ArepoMesh extent "));
		
		// preprocessing
		//ArepoMesh::ComputeVoronoiEdges();
		ArepoMesh::ComputeQuantityBounds();
		ArepoMesh::CalculateMidpoints();
		//ArepoMesh::LimitCellDensities();
		
		// copy individual allocation factors for auxiliary mesh
		AuxMesh.Indi.AllocFacNdp = AUXMESH_ALLOC_SIZE / 2;
		AuxMesh.Indi.AllocFacNdt = AUXMESH_ALLOC_SIZE;
		AuxMesh.Indi.AllocFacNvf = AUXMESH_ALLOC_SIZE;
		
		// TODO: temp units
		unitConversions[TF_VAL_DENS]   = All.UnitDensity_in_cgs / MSUN_PER_PC3_IN_CGS;
		unitConversions[TF_VAL_UTHERM] = All.UnitEnergy_in_cgs;
		
		IF_DEBUG(cout << "unitConv[dens]   = " << unitConversions[TF_VAL_DENS] << endl);
		IF_DEBUG(cout << "unitConv[utherm] = " << unitConversions[TF_VAL_UTHERM] << endl);
		
		// debugging
		//ArepoMesh::DumpMesh();
}

void ArepoMesh::LocateEntryCellBrute(const Ray &ray)
{
		// note: using the brute force search is O(N_rays * NumPart) - not good
		Point hitbox  = ray(ray.min_t);
		
		int DP_ID;
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
		if (DP_ID >= N_gas)
				cout << " LocateEntryCellBrute is [Local Ghost] DP_ID = " << DP_ID << " (N_gas = " << N_gas << ")" << endl;
#endif
		
		IF_DEBUG(cout << " brute DP_ID = " << DP_ID << " DP[DP_ID].index = " << DP[DP_ID].index 
									<< " dist = " << minDist << " (DP.x = " << DP[DP_ID].x << " DP.y = " 
									<< DP[DP_ID].y << " DP.z = " << DP[DP_ID].z << ")" << endl);
		
		ray.index = DP_ID;
		ray.task  = DP[DP_ID].task;
		
		// verify task assignment
		if (ray.task < 0 || ray.task >= NTask) {
				cout << "ERROR! ray has bad task=" << ray.task << endl;
				endrun(1115);
		}
		
}

// set ray.min_t and ray.max_t bounds first

void ArepoMesh::LocateEntryCell(const Ray &ray)
{
		Point hitbox  = ray(ray.min_t);
		Point exitbox = ray(ray.max_t);
		
		IF_DEBUG(cout << " ray starts at x = " << hitbox.x << " y = " << hitbox.y << " z = " << hitbox.z << endl);
		IF_DEBUG(cout << " ray ends at   x = " << exitbox.x << " y = " << exitbox.y << " z = " << exitbox.z << endl);		
		
		// use peanokey to find domain and task for ray
		if (ray.task == -1) {
				// TODO
				ray.task = 0;
		}
		// TODO: exchange
		
		// use tree to find nearest gas particle (local only)
		double mindist, mindist2;
		
#ifndef USE_AREPO_TREEFIND_FUNC
		int use_periodic = 0;
		float hsml = 10.0; //TODO this is really the "guess" which should be good (small) so that not too many nodes are opened, but unknown if this is effectively the "max" so that if this is too small no point will be found?
		const int dp_min = ArepoMesh::FindNearestGasParticle(hitbox, hsml, &mindist, -1, use_periodic);
#endif
		
#ifdef USE_AREPO_TREEFIND_FUNC
		endrun(6233);
		mindist = -1;
		double pos[3];
		int found_global;
		pos[0] = hitbox.x;
		pos[1] = hitbox.y;
		pos[2] = hitbox.z;
		const int dp_min = 0; // TODO 
		//= ngb_treefind_nearest_local(pos,-1,&mindist,0,&found_global);
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
						endrun(1139);
				}
		}
		
		if (count > 10) {
				cout << "WARNING: LocateEntryCell iterated [" << count << "] times." << endl;
				cout << " tree found dp_min=" << dp_min << " x = " << DP[dp_min].x << " y = " << 
								DP[dp_min].y << " z = " << DP[dp_min].z << endl;
				cout << " refine ended on dp_new=" << dp_new << " x = " << DP[dp_new].x << " y = " << 
								DP[dp_new].y << " z = " << DP[dp_new].z << endl;
				cout << " ray (hunting for, hitbox) at x = " << hitbox.x << " y = " << hitbox.y <<
								" z = " << hitbox.z << endl;
		}
		
		// if we did not finish in a primary cell, check that we iterated over at least one neighbor
		if (dp_new >= N_gas && !count) {
				cout << "ERROR: Refined entry tree search ended in ghost but count=0" << endl;
				endrun(1107);
		}
		
		ray.index = dp_new;
		ray.task  = DP[dp_new].task;
		
		// verify task assignment
		if (ray.task < 0 || ray.task >= NTask) {
				cout << "ERROR! ray has bad task=" << ray.task << endl;
				endrun(1115);
		}
		
}

void ArepoMesh::VerifyPointInCell(int dp, Point &pos)
{
		Vector celldist(pos.x - DP[dp].x, pos.y - DP[dp].y, pos.z - DP[dp].z);
		
		double dist2point = celldist.LengthSquared();
		
		// exactly on DP point, avoid divide by zero
		if (dist2point == 0.0)
				return;
				
		double dist2point_min = dist2point;
		
		// check neighbors for a closer tetra midpoint
		const int start_edge = midpoint_idx[dp].first;
		const int n_edges    = midpoint_idx[dp].second;
		
		int c = -1;
		
		for (int i=0; i < n_edges; i++)
		{
				const int dp_neighbor = opposite_points[start_edge + i];
				
				celldist = Vector(pos.x - DP[dp].x, pos.y - DP[dp].y, pos.z - DP[dp].z);
				
				dist2point = celldist.LengthSquared();
				
				if (dist2point < dist2point_min) {
						dist2point_min = dist2point;
						c = dp_neighbor;
				}
		}
	
		// check
		if (dist2point/dist2point_min > 1+INSIDE_EPS) {
				cout << "VerifyPointInCell FAILED! dp = " << dp << " pos.x = " << pos.x << " pos.y = " << pos.y
				     << " pos.z = " << pos.z << endl;
				endrun(1129);
		}
		
		IF_DEBUG(cout << "VerifyPointInCell Passed! dp = " << dp << " pos.x = " << pos.x << " pos.y = " << pos.y
				          << " pos.z = " << pos.z << endl);

}

int ArepoMesh::FindNearestGasParticle(Point &pt, float hsml, double *mindist, int guess, int use_periodic)
{
	// based on ngbtree_walk.c:ngb_treefind_variable() (no MPI)
	int node, nearest, p;
	struct NgbNODE *current;
	double dx, dy, dz, cur_mindist, xtmp, ytmp, ztmp;
	float search_min[3], search_max[3], search_max_Lsub[3], search_min_Ladd[3];
	
	// search bounds
	for(i = 0; i < 3; i++)
	{
		search_min[i] = searchcenter[i] - hsml;
		search_max[i] = searchcenter[i] + hsml;
	}

	search_max_Lsub[0] = search_max[0] - boxSize_X;
	search_max_Lsub[1] = search_max[1] - boxSize_Y;
	search_max_Lsub[2] = search_max[2] - boxSize_Z;

	search_min_Ladd[0] = search_min[0] + boxSize_X;
	search_min_Ladd[1] = search_min[1] + boxSize_Y;
	search_min_Ladd[2] = search_min[2] + boxSize_Z;

	// starting node
	node = Ngb_MaxPart;
	
	if (guess >= 0) {
		nearest = guess;
	} else {
		// pick random gas particle for guess of the min distance (why?)
		//nearest = floor(get_random_number(SelRnd++) * N_gas);
		nearest = (int)floor(N_gas/2.0);
	}
	
	if (use_periodic) {
		dx = NGB_PERIODIC_LONG_X(P[nearest].Pos[0] - pt.x);
		dy = NGB_PERIODIC_LONG_Y(P[nearest].Pos[1] - pt.y);
		dz = NGB_PERIODIC_LONG_Z(P[nearest].Pos[2] - pt.z);
	} else {
		dx = fabs(P[nearest].Pos[0] - pt.x);
		dy = fabs(P[nearest].Pos[1] - pt.y);
		dz = fabs(P[nearest].Pos[2] - pt.z);
	}
	cur_mindist = sqrt(dx * dx + dy * dy + dz * dz);

	while(node >= 0)
	{
		if(node < Ngb_MaxPart)  // single particle
		{
			p = node;
			node = Ngb_Nextnode[node];

			if(P[p].Type > 0) // not gas particle
				continue;

			/*
			if (use_periodic)
				dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - pt.x);
			else
				dx = fabs(P[p].Pos[0] - pt.x);
			if(dx > cur_mindist)
					continue;
					
			if (use_periodic)
				dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - pt.y);
			else
				dy = fabs(P[p].Pos[1] - pt.y);
			if(dy > cur_mindist)
					continue;
				
			if (use_periodic)
				dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - pt.z);
			else
				dz = fabs(P[p].Pos[2] - pt.z);
			if(dz > cur_mindist)
					continue;
			*/
			
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
			if(curdist2 > cur_mindist * cur_mindist)
				continue;

			cur_mindist = sqrt(curdist2);
			nearest = p;
		}
		else if(node < Ngb_MaxPart + Ngb_MaxNodes) // internal node
		{
			current = &Ngb_Nodes[node];

			// in case the node can be discarded
			node = current->u.d.sibling;

			// first quick tests along the axes (OLD)
			/*
			double test_dist = cur_mindist + 0.5 * current->len;
			
			if (use_periodic)
				dx = NGB_PERIODIC_LONG_X(current->center[0] - pt.x);
			else
				dx = fabs(current->center[0] - pt.x);
			if(dx > test_dist)
					continue;
					
			if (use_periodic)
				dy = NGB_PERIODIC_LONG_Y(current->center[1] - pt.y);
			else
				dy = fabs(current->center[1] - pt.y);
			if(dy > test_dist)
					continue;
					
			if (use_periodic)
				dz = NGB_PERIODIC_LONG_Z(current->center[2] - pt.z);
			else
				dz = fabs(current->center[2] - pt.z);
			if(dz > test_dist)
					continue;

			// now test against the minimal sphere enclosing everything
			test_dist += FACT1 * current->len;
			if(dx * dx + dy * dy + dz * dz > test_dist * test_dist)
					continue;
			*/
			
			// NEW
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
	}

	*mindist = cur_mindist;

	if (nearest < 0 || nearest > N_gas) {
		cout << "ERROR: FindNearestGasParticle nearest=" << nearest << " out of bounds." << endl;
		endrun(1118);
	}
	
	return nearest;
}

// get primary hydro ID - handle local ghosts
inline int ArepoMesh::getSphPID(int dp_id)
{
		int SphP_ID = -1;
		
		if (dp_id >= 0 && dp_id < N_gas)
				SphP_ID = dp_id;
		else if (dp_id >= N_gas)
				SphP_ID = dp_id - N_gas;
		
		if (SphP_ID < 0)
				endrun(1135);
}

bool ArepoMesh::AdvanceRayOneCellNew(const Ray &ray, float *t0, float *t1, 
																		 int previous_cell, Spectrum &Lv, Spectrum &Tr)
{
		double min_t = MAX_REAL_NUMBER;
		int qmin = -1; // next
		
		// verify task
		if (ray.task != ThisTask) {
				cout << "[" << ThisTask << "] ERROR! ray.task = " << ray.task << " differs from ThisTask!" << endl;
				endrun(1138);
		}
	
#ifdef DEBUG_VERIFY_INCELL_EACH_STEP
		// verify ray is where we expect it
		int dp = ray.index;
		Point pos = ray(ray.min_t);
		ArepoMesh::VerifyPointInCell(dp,pos);
#endif
		
		int SphP_ID = getSphPID(ray.index);

		// tetra position
		const Vector cellp(DP[ray.index].x,DP[ray.index].y,DP[ray.index].z);

		// just for verbose reporting and for the midp vector
		Point hitbox  = ray(*t0);
		Point exitbox = ray(*t1);
		
		IF_DEBUG(hitbox.print(" hitbox "));
		IF_DEBUG(exitbox.print(" exitbox "));
		IF_DEBUG(cellp.print(" cellp (DP.xyz) "));
		
		// alternative connectivity
		const pair<int,int> edge = midpoint_idx[ray.index];
		
		for (int i=edge.second-1; i >= 0; i--)
		{
				// skip face we arrived through, if any
				if (opposite_points[edge.first + i] == previous_cell && previous_cell != -1)
						continue;
				
				IF_DEBUG(cout << " checking face[" << i << "] midpoints.x = " << midpoints[edge.first + i].x
				              << " midpoints.y = " << midpoints[edge.first + i].y << " midpoints.z = " 
											<< midpoints[edge.first + i].z << " opposite_point = " 
											<< opposite_points[edge.first + i] << endl);
		
				// midpoint (c)
				const Vector midp(midpoints[edge.first + i].x - hitbox.x,
													midpoints[edge.first + i].y - hitbox.y,
													midpoints[edge.first + i].z - hitbox.z); //TODO: verify hitbox not hitcell
				
				// vector pointing to the outside, normal to a voronoi face of the cell (q)
				const Vector norm = midpoints[edge.first + i] - cellp;
		
				IF_DEBUG(midp.print(" midp (midpoints-hitbox) "));
				IF_DEBUG(norm.print(" norm (midpoints-cellp) "));
		
				// find intersection of ray with this face
				double dotprod1 = Dot( ray.d, norm ); //ddotq
				double dotprod2 = Dot( midp,  norm ); //cdotq
				
				// check if ray is aligned on face (e.g. backgroundgrid)
				if (dotprod1 == 0 && dotprod2 == 0)
						continue;
				
				if (dotprod1 > 0) {
						double t = dotprod2 / dotprod1; // infinite line/plane intersection test
						
						IF_DEBUG(cout << " i[" << i << "] dotprod>0 and t = " << t << " (min=" << (ray.min_t-*t0) << ")" << endl);
						
						if (t > (ray.min_t-*t0) && t < min_t) {
								min_t = t;
								qmin = opposite_points[edge.first + i];
								IF_DEBUG(cout << " intersection t = " << t << " setting new min_t, qmin (next DP) = " << qmin << endl);
						}
				} else {
						//cout << "Warning: Point on wrong side of face." << endl;
						// check if numerical error has put it on the wrong side of the face
						if (dotprod2 > 0) {
								//endrun(1137);
								
								//TODO: check and allow if we are really close to a face
						}
				}
		}
		
		// check if exiting box and failed to exit a face
		if (qmin == -1) {
				Point exitcell = ray(*t0 + min_t);
				
				if (!extent.Inside(exitcell)) { // && ray.index >= N_gas
						// set intersection with box face to allow for final contribution to ray
						IF_DEBUG(cout << " failed to intersect face, exitcell outside box, ok!" << endl);
						min_t = ray.max_t;
						
						// fake exit face
						qmin = 0;
				}
		}
		
		//TODO: temp disable adding any contribution from ghosts to rays
		bool addFlag = true;
		if (qmin != -1 && ray.index >= N_gas)
				addFlag = false;
		
		// check for proper exit point
		if (qmin != -1)
		{
				// clamp min_t to avoid integrating outside the box
				IF_DEBUG(cout << " min_t = " << min_t << " (t1=" << *t1 << " t0=" << *t0 << ")" << endl);
				min_t = Clamp(min_t,0.0,(*t1-*t0));
				
				if (addFlag)
				{ //TODO
				
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
						
						IF_DEBUG(midpt.print(" onestep midp rel sphcen "));
						
						double rho    = SphP[SphP_ID].Density;
						double utherm = SphP[SphP_ID].Utherm;
						
						// optical depth: treat as constant over entire cell
						Spectrum stepTau(0.0);
						
						// use gradients if requested
						if (Config.useDensGradients) {
								rho += Dot(sphDensGrad,midpt);
								stepTau += transferFunction->sigma_t() * rho * len;
						} else {
								// always use gradient for tau calculation
								stepTau += transferFunction->sigma_t() * (rho + Dot(sphDensGrad,midpt)) * len;
						}
						//if (Config.useUthermGradients) //TODO: gradient (MATERIALS)
						//		utherm += Dot(sphUthermGrad,midpt);
						
						// reduce transmittance for optical depth
						//Tr *= Exp(-stepTau); //TODO
						
						// TODO: sub-step length should be adaptive based on gradients
						// TODO: should only substep if the TF will evaluate nonzero somewhere inside

						// TODO: fix non-sub-stepping!!
						
						// if not sub-stepping then set default
						if (!viStepSize)
							viStepSize = len;
							
						// sub-stepping: variable 
						//int nSamples = (int)ceilf(len_sample / viStepSize);
						//double fracstep = 1.0 / nSamples;
						//double halfstep = 0.5 / nSamples;
						
						// sub-stepping: strict in world space
						Vector prev_sample_pt(ray(ray.depth * viStepSize));
						Vector norm_sample(exitcell - prev_sample_pt);

						int nSamples = (int)floorf(norm_sample.Length() / viStepSize);
						double fracstep = 1.0 / nSamples;
						norm_sample = Normalize(norm_sample);
						
						IF_DEBUG(prev_sample_pt.print(" prev_sample_pt "));
						
						IF_DEBUG(cout << " sub-stepping len = " << len << " nSamples = " << nSamples 
													<< " (step = " << len/nSamples << ")" << endl);
						
						for (int i = 0; i < nSamples; ++i) {
								//Vector midpt(hitcell[0] + (i*fracstep+halfstep) * norm[0],
								//						 hitcell[1] + (i*fracstep+halfstep) * norm[1],
								//						 hitcell[2] + (i*fracstep+halfstep) * norm[2]);
														 
								Vector midpt(prev_sample_pt[0] + ((i+1)*viStepSize) * norm_sample[0],
														 prev_sample_pt[1] + ((i+1)*viStepSize) * norm_sample[1],
														 prev_sample_pt[2] + ((i+1)*viStepSize) * norm_sample[2]);				 
								
								IF_DEBUG(midpt.print(" substep midpt "));
								
								double rho = nnInterpScalar(SphP_ID, ray.index, midpt);
								
								IF_DEBUG(cout << "  segment[" << i << "] fractrange [" << (i*fracstep) << "," 
															<< (i*fracstep)+fracstep << "] rho = " << SphP[SphP_ID].Density
															<< " rho nni = " << rho << endl);
								
								// pack cell quantities for TF
								float vals[TF_NUM_VALS];
								
								vals[TF_VAL_DENS]        = (float) rho;
								vals[TF_VAL_UTHERM]      = (float) SphP[SphP_ID].Utherm; //TODO: gradient (MATERIALS)
								vals[TF_VAL_PRES]        = (float) SphP[SphP_ID].Pressure; //TODO: gradient
								vals[TF_VAL_ENERGY]      = (float) SphP[SphP_ID].Energy;
								
								vals[TF_VAL_VEL_X]       = (float) P[SphP_ID].Vel[0]; //TODO: gradients
								vals[TF_VAL_VEL_Y]       = (float) P[SphP_ID].Vel[1];
								vals[TF_VAL_VEL_Z]       = (float) P[SphP_ID].Vel[2];
								//vals[TF_VAL_VEL_DIV]     = (float) SphP[SphP_ID].DivVel;
								//vals[TF_VAL_VEL_CURL]    = (float) SphP[SphP_ID].CurlVel;
								
								//vals[TF_VAL_POTENTIAL]   = (float) P[SphP_ID].Potential;
								
#ifdef METALS
								vals[TF_VAL_METALLICITY] = (float) SphP[SphP_ID].Metallicity;
#endif
#ifdef COOLING
								vals[TF_VAL_NE]          = (float) SphP[SphP_ID].Ne;
#endif
#ifdef USE_SFR
								vals[TF_VAL_SFR]         = (float) SphP[SphP_ID].Sfr;
#endif

								// apply TF to integrated (total) quantities (only appropriate for constant?)
								if (Config.projColDens) {
										rho    *= len;
										utherm *= len;
								}
						
								// compute emission-only source term using transfer function
								Lv += Tr * transferFunction->Lve(vals) /* * fracstep*/;
								
							// update previous sample point marker
							ray.depth++;
						}
		
				} //addFlag //TODO
				
				// update ray: transfer to next voronoi cell (possibly on different task)
				ray.task  = DP[qmin].task;
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
						cout << " ray pos: " << hitbox.x << " " << hitbox.y << " " << hitbox.z << endl;
						cout << " DP[ray.index] pos: " << DP[ray.index].x << " " << DP[ray.index].y 
								 << " " << DP[ray.index].z << endl << endl;
						endrun(1130);
				}
		}
		
		return true;		
}

// natural neighbor interpolation on scalar field at position pt inside Voronoi cell SphP_ID
inline double ArepoMesh::nnInterpScalar(int SphP_ID, int DP_ID, Vector &pt)
{
#ifdef NATURAL_NEIGHBOR_INTERP
		int tlast = 0;
		
		// check degenerate point in R3, immediate return
		if (fabs(pt.x - DP[DP_ID].x) <= INSIDE_EPS &&
				fabs(pt.y - DP[DP_ID].y) <= INSIDE_EPS &&
				fabs(pt.z - DP[DP_ID].z) <= INSIDE_EPS)
				return SphP[SphP_ID].Density;
				
		// list of neighbors
		const int start_edge = midpoint_idx[DP_ID].first;
		const int n_edges    = midpoint_idx[DP_ID].second;
		
		int *dp_neighbors   = new int[n_edges+1];
		int *sphp_neighbors = new int[n_edges+1];
		
		for (int i=0; i < n_edges; i++) {
				dp_neighbors[i]   = opposite_points[start_edge + i];
				sphp_neighbors[i] = getSphPID(opposite_points[start_edge + i]);
		}
		// add parent to list
		dp_neighbors[n_edges]   = DP_ID;
		sphp_neighbors[n_edges] = SphP_ID;
		
#ifdef DEBUG
		cout << " dp_neighbors: ";
		for (int i=0; i < n_edges; i++)
				cout << " " << opposite_points[start_edge + i] << "(sphp=" << 
								getSphPID(opposite_points[start_edge + i]) << ")";
		cout << " " << DP_ID << " (sphp=" << SphP_ID << ")";
		cout << " [" << n_edges+1 << " neighbors]" << endl;
#endif
		
		//P[i] = cell to be split ("active") or in our case the parent of the interpolate pt
		//j = N_gas + count; ("new P/SphP id")
		//P[j] and SphP[j] are new, i is old
		//jj = ref_SphP_dp_index[i]; //dp index of SphP point i
		int jj = DP_ID;
		
		// clear auxiliary mesh (keep super-tetra)
		/*
		IF_DEBUG(cout << "DeRefMesh.Nvf = " << DeRefMesh.Nvf << " Ndp = " << DeRefMesh.Ndp << " Ndt = " <<
		         DeRefMesh.Ndt << endl);
		if (DeRefMesh.Ndp > 5) {
			memset(DeRefMesh.VF,0,DeRefMesh.MaxNvf * sizeof(face));
			memset(DeRefMesh.DP,0,(DeRefMesh.MaxNdp - 5) * sizeof(point));
			memset(DeRefMesh.DT,0,DeRefMesh.MaxNdt * sizeof(tetra));
			memset(DeRefMesh.DTC,0,DeRefMesh.MaxNdt * sizeof(tetra_center));
			
			for(int k = 0; k < DeRefMesh.MaxNdt; k++)
			  DeRefMesh.DTF[k] = 0;
				
			// set aux mesh sizes to zero
			DeRefMesh.Ndp = 5;
			DeRefMesh.Ndt = 0;
			DeRefMesh.Nvf = 0;
		} */
		
		// recreate the voronoi cell of the parent's neighbors (auxiliary mesh approach)
		initialize_and_create_first_tetra(&DeRefMesh);
		
		DeRefMesh.DTC = static_cast<tetra_center*>(
										mymalloc_movable(&DeRefMesh.DTC, "AuxDTC", DeRefMesh.MaxNdt * sizeof(tetra_center)));
		DeRefMesh.DTF = static_cast<char*>(
										mymalloc_movable(&DeRefMesh.DTF, "AuxDTF", DeRefMesh.MaxNdt * sizeof(char)));
		
		// construct new auxiliary mesh around pt
		for(int k = 0; k < n_edges+1; k++)
		{
				int q = dp_neighbors[k];
		
				if(DeRefMesh.Ndp + 2 >= DeRefMesh.MaxNdp)
				{
						DeRefMesh.Indi.AllocFacNdp *= ALLOC_INCREASE_FACTOR;
						DeRefMesh.MaxNdp = (int)DeRefMesh.Indi.AllocFacNdp;
						DeRefMesh.DP -= 5;
						DeRefMesh.DP = static_cast<point*>
													 (myrealloc_movable(DeRefMesh.DP, (DeRefMesh.MaxNdp + 5) * sizeof(point)));
						DeRefMesh.DP += 5;
						//endrun(1157);
				}
				
				DeRefMesh.DP[DeRefMesh.Ndp] = Mesh.DP[q];
		    set_integers_for_point(&DeRefMesh, DeRefMesh.Ndp);
		    tlast = insert_point(&DeRefMesh, DeRefMesh.Ndp, tlast);
		    DeRefMesh.Ndp++;
				
				//IF_DEBUG(cout << " inserted neighbor k=" << k << " new tlast=" << tlast << " totnum=" << 
				//				 DeRefMesh.Ndp << endl);
		}
		
		// compute old circumcircles and volumes
		compute_circumcircles(&DeRefMesh);
		
		double *dp_old_vol = new double[DeRefMesh.Ndp+1];
		double *dp_new_vol = new double[DeRefMesh.Ndp+1];
		
		//computeAuxVolumes(dp_old_vol);
		derefine_refine_compute_volumes(dp_old_vol);
		
#ifdef DEBUG
		cout << " old volumes:";
		for (int k = 0; k < DeRefMesh.Ndp; k++) {
				cout << " [k=" << k << "] " << dp_old_vol[k];
		}
		cout << endl;
#endif		
		
		// add interpolate point
		DeRefMesh.DP[DeRefMesh.Ndp].x = pt.x;
		DeRefMesh.DP[DeRefMesh.Ndp].y = pt.y;
		DeRefMesh.DP[DeRefMesh.Ndp].z = pt.z;
		DeRefMesh.DP[DeRefMesh.Ndp].ID = -1; //unused
		set_integers_for_point(&DeRefMesh, DeRefMesh.Ndp);
		tlast = insert_point(&DeRefMesh, DeRefMesh.Ndp, tlast);
		DeRefMesh.Ndp++;
		
		//IF_DEBUG(cout << " inserted interpolate new tlast=" << tlast << " totnum=" << DeRefMesh.Ndp << endl);
		
		// compute new circumcircles and volumes
		compute_circumcircles(&DeRefMesh);
		//computeAuxVolumes(dp_new_vol);
		derefine_refine_compute_volumes(dp_new_vol);
		
#ifdef DEBUG
		cout << " new volumes:";
		for (int k = 0; k < DeRefMesh.Ndp; k++) {
				cout << " [k=" << k << "] " << dp_new_vol[k];
		}
		cout << endl;
#endif

		// calculate scalar value based on neighbor values and area fraction weights
		double *weights = new double[DeRefMesh.Ndp-1];
		double invInterpVol = 1.0 / dp_new_vol[DeRefMesh.Ndp-1];
		
		for (int k=0; k < DeRefMesh.Ndp-1; k++)
			weights[k] = (dp_old_vol[k]-dp_new_vol[k]) * invInterpVol;
			
#ifdef DEBUG
		cout << " weights (vals):";
		for (int k=0; k < DeRefMesh.Ndp-1; k++)
				cout << " [k=" << k << "]" << weights[k] << " (" << sphp_neighbors[k] << ")";
		cout << endl;
#endif

		double val = 0.0;
		for (int k=0; k < DeRefMesh.Ndp-1; k++) {
				val += SphP[sphp_neighbors[k]].Density * weights[k];
		}
		
		// free
		delete weights;
		delete dp_new_vol;
		delete dp_old_vol;
		delete sphp_neighbors;
		delete dp_neighbors;
		
		// free aux mesh
		myfree(DeRefMesh.DTF);
		myfree(DeRefMesh.DTC);
		DeRefMesh.DTC = NULL;
		myfree(DeRefMesh.DT);
		myfree(DeRefMesh.DP - 5);
		myfree(DeRefMesh.VF);
#endif //NNI

#ifndef NATURAL_NEIGHBOR_INTERP
		// old:
		Vector sphCen(     SphP[SphP_ID].Center[0],   SphP[SphP_ID].Center[1],   SphP[SphP_ID].Center[2]);
		Vector sphDensGrad(SphP[SphP_ID].Grad.drho[0],SphP[SphP_ID].Grad.drho[1],SphP[SphP_ID].Grad.drho[2]);
    pt -= sphCen;	// make relative to cell center	
		double val   = SphP[SphP_ID].Density + Dot(sphDensGrad,pt);
		
		//IF_DEBUG(cout << "NNI: " << val << " LinearGrad: " << val_old << " frac diff = " << 
		//				fabs(val-val_old)/val << endl);
#endif

		return val;
}

inline void ArepoMesh::computeAuxVolumes(double *vol)
{
		int i, bit, nr;

		for(i = 0; i < AuxMesh.Ndp; i++)
				vol[i] = 0;

		Edge_visited = static_cast<unsigned char*>(
									mymalloc_movable(&Edge_visited, "Edge_visited", AuxMesh.Ndt * sizeof(unsigned char)));

		for(i = 0; i < AuxMesh.Ndt; i++)
				Edge_visited[i] = 0;

		for(i = 0; i < AuxMesh.Ndt; i++)
		{
				if(AuxMesh.DT[i].t[0] < 0)	/* deleted ? */
						continue;

				bit = 1;
				nr  = 0;

				while(Edge_visited[i] != EDGE_ALL)
				{
						if((Edge_visited[i] & bit) == 0)
						derefine_refine_process_edge(&AuxMesh, vol, i, nr);

						bit <<= 1;
						nr++;
				}
		}

		myfree(Edge_visited);
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
				
						int SphP_ID = -1;
						
						if (DP[dp].index >= 0 && DP[dp].index < N_gas)
								SphP_ID = DP[dp].index;
						else if (dp >= N_gas)
								SphP_ID = DP[dp].index - N_gas;
								
						// valid cell?
						if (DP[dp].index < N_gas && SphP_ID >= 0) {
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
				int SphP_ID = -1;
				
				if (DP[i].index >= 0 && DP[i].index < N_gas)
						SphP_ID = DP[i].index;
				else if (i >= N_gas)
						SphP_ID = DP[i].index - N_gas;
				
				if (SphP_ID < 0)
						endrun(1133);
				
				mm.insert(make_pair(SphP_ID, i));
		}
		
		// set up mapping of DP id -> DP primary id
		for (int i=0; i < Ndp; i++)
		{
				int SphP_ID = -1;
				
				if (DP[i].index >= 0 && DP[i].index < N_gas)
						SphP_ID = DP[i].index;
				else if (i >= N_gas)
						SphP_ID = DP[i].index - N_gas;
		
				if (SphP_ID < 0)
						endrun(1134);
		
				// cell has no hydro quantities -> map to -1
				if (SphP_ID < 0) {
						IF_DEBUG(cout << "WARNING: CM i=" << i << " SphP_ID (neg) = " << SphP_ID << endl);
						primary_cells.push_back(-1);
				}
				
				// loop over all DP indices that share this SphP cell
				pair<mmi,mmi> dp_indices(mm.equal_range(SphP_ID));
				
				if (dp_indices.second == dp_indices.first)
						endrun(1131);
						
				// search for primary cell and record to map
				while (dp_indices.first != dp_indices.second)
				{
						if (dp_indices.first->second >= 0 && dp_indices.first->second < N_gas) {
								primary_cells.push_back(dp_indices.first->second);
								break;
						}
						dp_indices.first++;
				}
		}
		
		// verify size
		if (primary_cells.size() != Ndp)
				endrun(1132);
				
		// use VF array to generate connnections (reorganize it to index by point)
		multimap<int,int> conn;
		
		for (int i=0; i < Nvf; i++) {
				conn.insert(make_pair(VF[i].p1,VF[i].p2));
				conn.insert(make_pair(VF[i].p2,VF[i].p1));
		}
		
		// associate connections with cells
		for (int i=0; i < Ndp; i++)
		{
				int SphP_ID = -1;
				
				if (DP[i].index >= 0 && DP[i].index < N_gas)
						SphP_ID = DP[i].index;
				else if (i >= N_gas)
						SphP_ID = DP[i].index - N_gas;
						
				if (SphP_ID < 0)
						continue;
						
				const Vector cellp(DP[i].x,DP[i].y,DP[i].z);
				
				// find connections for this cell
				pair<mmi,mmi> dp_neighbors(conn.equal_range(i));
				
				if (dp_neighbors.first == dp_neighbors.second)
						endrun(1133);
						
				for (; dp_neighbors.first != dp_neighbors.second; dp_neighbors.first++)
				{
						const int dp_neighbor = dp_neighbors.first->second;
				
						int SphP_ID_n = -1;
						
						if (DP[dp_neighbor].index >= 0 && DP[dp_neighbor].index < N_gas)
								SphP_ID_n = DP[dp_neighbor].index;
						else if (dp_neighbor >= N_gas) {
								//cout << " N_gas=" << N_gas << " dp_neighbor=" << dp_neighbor << " DP[dp_neighbor].index="
								//     << DP[dp_neighbor].index << endl;
								SphP_ID_n = DP[dp_neighbor].index - N_gas;
						}
								
						//cout << "pri=" << setw(2) << i << " dp_neighbor=" << dp_neighbor << " SphID=" << SphP_ID_n << endl;
								
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
		
		for (int i=0; i < N_gas; i++) {
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
		valBounds[TF_VAL_DENS*3 + 2] = pmean / N_gas;
		
		valBounds[TF_VAL_UTHERM*3 + 0] = umin;
		valBounds[TF_VAL_UTHERM*3 + 1] = umax;
		valBounds[TF_VAL_UTHERM*3 + 2] = umean / N_gas;
		
		IF_DEBUG(cout << " Density min = " << valBounds[TF_VAL_DENS*3 + 0] 
									<< " max = " << valBounds[TF_VAL_DENS*3 + 1] 
									<< " mean = " << valBounds[TF_VAL_DENS*3 + 2] << endl);
									
		IF_DEBUG(cout << " Utherm  min = " << valBounds[TF_VAL_UTHERM*3 + 0] 
									<< " max = " << valBounds[TF_VAL_UTHERM*3 + 1] 
									<< " mean = " << valBounds[TF_VAL_UTHERM*3 + 2] << endl);
									
}

void ArepoMesh::ComputeVoronoiEdges()
{
		IF_DEBUG(cout << "ArepoMesh::ComputeVoronoiEdges()" << endl);
		
		// geometric conventions (voronoi.h)
		const int edge_start[6]     = { 0, 0, 0, 1, 1, 2 };
		const int edge_end[6]       = { 1, 2, 3, 2, 3, 3 };
		const int edge_opposite[6]  = { 3, 1, 2, 3, 0, 1 };
		const int edge_nexttetra[6] = { 2, 3, 1, 0, 2, 0 };
		
		int Nel    = 0;
		int MaxNel = 10 * Nvf;
		int startingTetra = -1;
		
		tetra *prev, *next;
		int i_face,i,j,k,l,m,ii,jj,kk,ll,nn,pp,tt,nr,nr_next;
		
		int dp1,dp2;		
		
		// allocate
		EdgeList     = new int[MaxNel];
		Nedges       = new int[Nvf];
		NedgesOffset = new int[Nvf];
		Edge_visited = new unsigned char[T->Ndt];
		
		// zero
		for(i = 0; i < Nvf; i++) {
				Nedges[i]     = 0;
    }

		// for each hydro point, find a tetra with it as a vertex (rough)
		//for(i = 0; i < T->Ndt; i++) {
		//		for(j = 0; j < DIMS + 1; j++) {
		//				if(DP[DT[i].p[j]].index >= 0 && DP[DT[i].p[j]].index < N_gas) // local, non-ghost, TODO: removed task
		//						whichtetra[DP[DT[i].p[j]].index] = i;
		//		}
    //}
		
		for (i_face = 0; i_face < Nvf; i_face++) {

				// locate tetra containing both points
				dp1 = VF[i_face].p1;
				dp2 = VF[i_face].p2;
				
				// skip voronoi faces that have vertices in the boundary tetra
				if (dp1 == -5 || dp2 == -5) {
						IF_DEBUG(cout << " VF[" << i_face << "] skipping! dp1 = " << dp1 << " dp2 = " << dp2 << endl << endl);
						continue;
				}
		
				// TODO: do this some other way, right now its Order(Nvf*Ndt) - not good
				for (i = 0; i < T->Ndt; i++) {
						if (DT[i].p[0] != dp1 && DT[i].p[1] != dp1 && DT[i].p[2] != dp1 && DT[i].p[3] != dp1)
						    continue;
						if (DT[i].p[0] != dp2 && DT[i].p[1] != dp2 && DT[i].p[2] != dp2 && DT[i].p[3] != dp2)
								continue;
						if (DT[i].t[0] >= 0) // deleted?
								startingTetra = i;
				}
				
				if (startingTetra < 0) {
						cout << "ComputeVoronoiEdges(" << i_face << ") ERROR! startingTetra = " << startingTetra << endl;
						endrun(1105);
				}
				
				// set starting tetra
				tetra *t = &DT[startingTetra];
				
				if (t->t[0] < 0) {
						cout << "ComputeVoronoiEdges(" << i_face << ") ERROR! Starting t->t[0] = " << t->t[0] << endl;
						endrun(1104);				
				}
				
				IF_DEBUG(cout << " VF[" << i_face << "] startingTetra = " << startingTetra << " (" << (t-DT) << ")" << endl);
		
				// set edge number to traverse around
				for (i = 0; i < 4; i++) {
						if (DT[startingTetra].p[i] == dp1) ii = i;
						if (DT[startingTetra].p[i] == dp2) jj = i;
				}
				
				if (ii > jj) swap(ii,jj);
				
				if (ii==0 && jj==1) nr = 0;
				if (ii==0 && jj==2) nr = 1;
				if (ii==0 && jj==3) nr = 2;
				if (ii==1 && jj==2) nr = 3;
				if (ii==1 && jj==3) nr = 4;
				if (ii==2 && jj==3) nr = 5;
				
				IF_DEBUG(cout << " starting ii = " << ii << " jj = " << jj << " nr = " << nr << endl);
				
				// zero
				for(i = 0; i < T->Ndt; i++)
						Edge_visited[i] = 0;
				
				i = edge_start[nr];
				j = edge_end[nr];
				k = edge_opposite[nr];
				l = edge_nexttetra[nr];

				Edge_visited[startingTetra] |= (1 << nr);
				
				pp = startingTetra;
				prev = t;

				do
				{
						nn   = prev->t[l];
						next = &DT[nn];
						
						Nedges[i_face]++;
						
						if(Nel >= MaxNel) {
								cout << "ERROR: Nel > MaxNel" << endl;
								endrun(1100);
						}
						
						EdgeList[Nel] = (int)(prev - DT); //i.e. element offset in DT/DTC [0,Ndt-1]
						IF_DEBUG(cout << " adding i = " << i_face << " EdgeList[" << Nel << "] = " << EdgeList[Nel] 
								          << " (num=" << Nedges[i_face] << ")" << endl);
						Nel++;
						
						// verify 3/4 (k,i,j) of tetra vertices are shared with the next tetra (opposite p[l])
						for(m = 0, ll = ii = jj = -1; m < 4; m++)
						{
								if(next->p[m] == prev->p[k])
									ll = m;
								if(next->p[m] == prev->p[i])
									ii = m;
								if(next->p[m] == prev->p[j])
									jj = m;
						}

						if(ll < 0 || ii < 0 || jj < 0) {
								cout << "ERROR: Inconsistency. ll = " << ll << " ii = " << ii << " jj = " << jj << " kk = " << kk
								     << " l = " << l << " i = " << i << " j = " << j << " k = " << k << endl;
								cout << "  prev p[0] = " << prev->p[0] << " p[1] = " << prev->p[1] << " p[2] = " 
								     << prev->p[2] << " p[3] = " << prev->p[3] << endl;
								cout << "  prev t[0] = " << prev->t[0] << " t[1] = " << prev->t[1] << " t[2] = " 
								     << prev->t[2] << " t[3] = " << prev->t[3] << endl;						 
								cout << "  next p[0] = " << next->p[0] << " p[1] = " << next->p[1] << " p[2] = " 
								     << next->p[2] << " p[3] = " << next->p[3] << endl;
								endrun(1101);
						}

						kk = 6 - (ll + ii + jj);
						
						// flag visited edge with bitmask
						for(nr_next = 0; nr_next < 6; nr_next++) {
								if((edge_start[nr_next] == ii && edge_end[nr_next] == jj) ||
									 (edge_start[nr_next] == jj && edge_end[nr_next] == ii))
								{
										if((Edge_visited[nn] & (1 << nr_next)) && next != t) {
												cout << "ERROR: Inconsistency (2)." << endl;
												endrun(1102);
										}

										//bitwise OR itself with 1*nr_next^2 (see EDGE_X)
										Edge_visited[nn] |= (1 << nr_next);
										break;
								}
						}
						
						IF_DEBUG(cout << " moving to tetra = " << nn << endl);
						
						prev = next;
						pp = nn;
						i = ii;
						l = ll;
						j = jj;
						k = kk;

				} while(next != t);

				IF_DEBUG(cout << " finished voronoi face i_face = " << i_face << " with Nedges = " 
											<< Nedges[i_face] << endl << endl);
				
				if (Nedges[i_face] < DIMS || Nedges[i_face] > 20) {
						cout << "ComputeVoronoiEdges(" << i_face << ") ERROR! Strange Nedges = " << Nedges[i_face] << endl;
						endrun(1103);
				}
		} //i_face
		
		/* // ugh 2D algorithm never going to work
		for(i = 0; i < N_gas; i++) {
				if(whichtetra[i] < 0)
						continue;

				qstart = q = &DT[whichtetra[i]];
				cout << " starting in tetra = " << whichtetra[i] << endl;
				do
				{
				
						Nedges[i]++;

						if(Nel >= MaxNel) {
								cout << "ERROR: Nel > MaxNel" << endl;
								return;
						}

						EdgeList[Nel] = q - DT; //i.e. element offset in DT/DTC [0,Ndt-1]
						cout << " adding i = " << i << " EdgeList[" << Nel << "] = " << EdgeList[Nel] << endl;
						Nel++;
						
						for(j = 0; j < DIMS+1; j++) {
								if(DP[q->p[j]].index == i) { // TODO: removed task
										cout << " found hydro point in tetra at index j = " << j << endl;
										break;
								}
						}

						// traverse
						k = j+1;

						if(k >= DIMS+1)
								k -= (DIMS+1);
						
						// move to adjacent tetra, opposite from p[k]
						cout << " moving to tetra = " << q->t[k] << endl;
						q = &DT[q->t[k]];
						
				} while(q != qstart);
				
				cout << " finished hydro point i = " << i << " with Nedges = " << Nedges[i] << endl << endl;
    }
		*/

		for(i = 1, NedgesOffset[0] = 0; i < Nvf; i++) {
				NedgesOffset[i] = NedgesOffset[i - 1] + Nedges[i - 1];
		}

}

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
		
		cout << endl << "SphP Hydro [" << N_gas << "]:" << endl;
		for (int i=0; i < N_gas; i++) {
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
		for (int i=0; i < N_gas; i++) {
				dc   = SphP[i].first_connection;
				cout << " SphP[" << setw(3) << i << "] DC.first = " << setw(2) << SphP[i].first_connection;
				
				do {
						next = DC[dc].next;
						cout << " " << " next = " << setw(2) << next;
						dc = next;
				} while (next != SphP[i].last_connection);
				
				cout << " DC.last = " << SphP[i].last_connection << endl;		
		}
		
		if (Nedges) {
			cout << endl << "Voronoi Edges (N_gas=" << N_gas << "):" << endl;
			for (int i=0; i < N_gas; i++) {
					cout << setw(3) << i << " Nedges = " << setw(2) << Nedges[i] << " NedgesOffset = " 
							 << setw(3) << NedgesOffset[i] << endl;
			}
		}
		
		cout << endl << "Primary_Cells and Midpoint_Idx (size=" << primary_cells.size() << "):" << endl;
		for (int i=0; i < primary_cells.size(); i++) {
				cout << "[" << setw(2) << i << "] primary id = " << primary_cells[i] 
				     << " edges start " << midpoint_idx[i].first << " num edges = " 
						 << midpoint_idx[i].second << endl;
		}
		
		cout << endl << "Midpoints and Opposite_Points (size=" << midpoints.size() << "):" << endl;
		for (int i=0; i < opposite_points.size(); i++) {
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
		IF_DEBUG(cout << "VoronoiEdges(" << i_face << ") Nedges=" << Nedges[i_face]
									<< " NedgesOffset = " << NedgesOffset[i_face] << endl);
		
		if (Nedges[i_face] <= 0 || Nedges[i_face] < DIMS || i_face < 0 || i_face >= Nvf) {
				IF_DEBUG(cout << "WARNING: Nedges[" << i_face << "] empty, degenerate, or out of bounds." << endl);
				return false;
		}
		
		int s_ind = EdgeList[NedgesOffset[i_face]];
		int n_ind;
		
		Point prev = Point(DTC[s_ind].cx, DTC[s_ind].cy, DTC[s_ind].cz);
		Point next;

		for (int i=1; i < Nedges[i_face]+1; i++) {				
				n_ind = EdgeList[(NedgesOffset[i_face] + i) % Nedges[i_face]];
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
