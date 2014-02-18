/*
 * arepo.cpp
 * dnelson
 */
 
#include "transform.h"
#include "spectrum.h"
#include "volume.h"
#include "transfer.h"

#include "arepo.h"
#include "util.h" // for numberOfCores()
#include "snapio.h"

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

		determine_compute_nodes();

		int num_threads = 1;
#if (NUM_THREADS > 1)
		//calling omp_set_num_threads overrides OMP_NUM_THREADS env variable
		omp_set_num_threads(NUM_THREADS);
		omp_set_dynamic(0);
		num_threads = omp_get_max_threads();
#endif
		
		cout << "AREPO ENABLED. (OpenMP NThreads = " << num_threads << ", MPI NTask = " << NTask 
		     << " ThisTask = " << ThisTask << ")" << endl;
}

void Arepo::Cleanup()
{
		MPI_Finalize();
		close_logfiles();
}

bool Arepo::LoadSnapshot()
{
    IF_DEBUG(cout << "Arepo::LoadSnapshot(" << snapFilename << ")." << endl);
		
		if( !Config.verbose )
			freopen("/dev/null","w",stdout); //hide arepo stdout	
		
		// set startup options
		WriteMiscFiles = 0;
		RestartSnapNum = -1;
		RestartFlag    = SUNRISE_CODE;

		strcpy(ParameterFile,paramFilename.c_str());

		// call arepo: run setup (also sets units)
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
		
		if( Config.totNumJobs >= 1 ) {
			// custom selective load snapshot, make/use maskfile for job splitting
			ArepoSnapshot arepoSnap( snapFilename );
			arepoSnap.read_ic();
		} else {		
			// load snapshot (GAS ONLY) with Arepo function
			read_ic(snapFilename.c_str(), 0x01);
		}
		
		if( Config.nTreeNGB )
		{
			// sampling based on tree searches only, custom (minimal) init
			int i, j;
			
			All.TotNumOfForces = 0;
			All.TopNodeAllocFactor = 0.08;
			All.TreeAllocFactor = 0.7;
			All.NgbTreeAllocFactor = 0.7;
			
			if(All.ComovingIntegrationOn)	//  change to new velocity variable
			{
				for(i = 0; i < NumPart; i++)
				{
					for(j = 0; j < 3; j++)
						P[i].Vel[j] *= sqrt(All.Time) * All.Time;	// for dm/gas particles, p = a^2 xdot
				}
			}
			
			domain_Decomposition();
			set_softenings();
			
			// will build tree
			Ngb_MaxPart = All.MaxPartSph;
			Ngb_MaxNodes = (int) (All.NgbTreeAllocFactor * All.MaxPartSph) + NTopnodes;
			ngb_treeallocate();
			ngb_treebuild(NumGas);
			
			//setup_smoothinglengths(); // build grav tree
			//cell density not set (by mesh), would need to include in snapshot file
			//perhaps most (all) fields not accessed, chance to really kill P/SphP memory usage
			//update_primitive_variables();	// to get pressure
			printf("Arepo tree loaded, returning control.\n");
		}
		else
		{
			// Voronoi mesh based sampling: call full init
			if (init() != SUNRISE_CODE) {
				cout << "Arepo::LoadSnapshot() ERROR: Arepo did not return successfully." << endl;
				return false;
			}
		}
		
		if( !Config.verbose ) { // restore stdout
			//TODO: switch between these automatically
			freopen("/dev/tty","w",stdout);
			
			//string fn = Config.imageFile + string(".out.txt");
			//freopen(fn.c_str(),"a",stdout); //return stdout to a file (SLURM)
		}		

		// illustris fof0 log density/temp
		/*
		cout << "Converting DENSITY TO LOG!" << endl;
		cout << "Converting TEMPERATURE TO LOG!" << endl;
		for(int i = 0; i < NumGas; i++)
		{
			SphP[i].Density = log10( SphP[i].Density );
			SphP[i].Utherm = log10( SphP[i].Utherm );
		} */

		// Data Processing (temporary spot)
		float gamma = 5.0 / 3.0;
		
		for( int i=0; i < NumGas; i++ )
		{
			// calculate entropy (leave in code units) (Pressure/prim vars set in init())
			SphP[i].OldMass = SphP[i].Pressure * All.cf_a3inv / pow(SphP[i].Density * All.cf_a3inv, gamma);
			
			if( SphP[i].OldMass < 0.0 ) {
				cout << "Error: Calculated [" << i << "] entropy = " << SphP[i].OldMass << endl;
				exit(1206);
			}
		}
		
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
		
		if( !All.BoxSize ) {
			cout << "Error: All.BoxSize=0, likely structure mismatch." << endl;
			exit(1179);
		}
		
		IF_DEBUG(extent.print(" ArepoMesh extent "));

#ifdef SPECIAL_BOUNDARY
		// override density of ID=-2 (inner spoon boundary cells) with density much higher than surrounding gas
		// put a TF near this value but below to highlight spoon boundary
		for(int i = 0; i < NumGas; i++)
		{
			if( P[i].ID == -2 )
				SphP[i].Density = 1000.0;
		}
#endif		
		
		// preprocessing
		ArepoMesh::ComputeQuantityBounds();
		//ArepoMesh::LimitCellDensities();
		
		ArepoMesh::setupAuxMeshes();
		ArepoMesh::precomputeTetraGrads();
		
		// TODO: temp units
		unitConversions[TF_VAL_DENS] = All.UnitDensity_in_cgs / MSUN_PER_PC3_IN_CGS;
		unitConversions[TF_VAL_TEMP] = All.UnitEnergy_in_cgs;
		
		IF_DEBUG(cout << "unitConv[dens]   = " << unitConversions[TF_VAL_DENS] << endl);
		IF_DEBUG(cout << "unitConv[utherm] = " << unitConversions[TF_VAL_TEMP] << endl);		
		
}

ArepoMesh::~ArepoMesh()
{				
#ifdef NATURAL_NEIGHBOR_INTERP
		// free aux meshes
		int numMeshes = numberOfCores();
		
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
}

void ArepoMesh::LocateEntryCellBrute(const Ray &ray)
{
		// note: using the brute force search is O(N_rays * NumPart) - not good
		Point hitbox  = ray(ray.min_t);
		
		int sphInd = -1;
		double minDist = MAX_REAL_NUMBER;
		Vector delta;
		double dist;
		
		for (int i=0; i < NumGas; i++) {
				delta.x = hitbox.x - P[i].Pos[0];
				delta.y = hitbox.y - P[i].Pos[1];
				delta.z = hitbox.z - P[i].Pos[2];
				
				dist = delta.PeriodicLengthSquared();
				
				if (dist < minDist) {
						minDist = dist;
						sphInd = i;
				}
		}
		
		IF_DEBUG(cout << " brute sphInd = " << sphInd
									<< " dist = " << minDist << " (P.x = " << P[sphInd].Pos[0] << " P.y = " 
									<< P[sphInd].Pos[1] << " P.z = " << P[sphInd].Pos[2] << ")" << endl);
		
		ray.index = sphInd;
		ray.task  = 0;
		
		// verify task assignment
		if (ray.task < 0 || ray.task >= NTask) {
				cout << "ERROR! ray has bad task=" << ray.task << endl;
				terminate("1115");
		}
		
}

// set ray.min_t and ray.max_t bounds first

void ArepoMesh::LocateEntryCell(const Ray &ray, int *prevEntryCell)
{
		double mindist;
	
		Point hitbox  = ray(ray.min_t);
		
#ifdef DEBUG
		Point exitbox = ray(ray.max_t);
		
		cout << " ray starts at x = " << hitbox.x << " y = " << hitbox.y << " z = " << hitbox.z << endl;
		cout << " ray ends at   x = " << exitbox.x << " y = " << exitbox.y << " z = " << exitbox.z << endl;
#endif		
		
		// TODO: use peanokey to find domain and task for ray
		if (ray.task == -1) {
				ray.task = 0;
		}
		// TODO: exchange
		
		// use tree to find nearest gas particle (local only)
		int dp_min = ArepoMesh::FindNearestGasParticle(hitbox, *prevEntryCell, &mindist);
		*prevEntryCell = dp_min;
		
		ray.index = dp_min;
		ray.task  = 0;
		
		// verify task assignment
		if (ray.task < 0 || ray.task >= NTask) {
				cout << "ERROR! ray has bad task=" << ray.task << endl;
				terminate("1115");
		}
		
}

void ArepoMesh::VerifyPointInCell(int parInd, Point &pos)
{		
		Vector celldist(pos.x - P[parInd].Pos[0], pos.y - P[parInd].Pos[1], pos.z - P[parInd].Pos[2]);
		double dist2point = celldist.PeriodicLengthSquared();
		
		if ( dist2point >= INSIDE_EPS )
		{
			// check all primary cells using periodic distance
			for (int i=0; i < NumGas; i++)
			{		
				celldist = Vector(pos.x - P[i].Pos[0], pos.y - P[i].Pos[1], pos.z - P[i].Pos[2]);
			
				if(celldist.PeriodicLengthSquared()/dist2point < 1 - INSIDE_EPS)
				{
					cout << "VerifyPointInCell FAILED! pt.x = " << setprecision(10) <<  pos.x << " pt.y = " << pos.y << " pt.z = " << pos.z << endl 
				       << " sphInd_cur = " << setw(3) << parInd << " P.x = " << P[parInd].Pos[0] << " P.y = " << P[parInd].Pos[1] 
				  		 << " P.z = " << P[parInd].Pos[2] << " (dist2point = " << dist2point << ")" << endl;
					cout << " sphInd_min = " << setw(3) << i << " P.x = " << P[i].Pos[0] << " P.y = " << P[i].Pos[1] 
				       << " P.z = " << P[i].Pos[2] << " (dist2point = " << celldist.PeriodicLengthSquared() << ")" << endl;
					terminate("1129");
				}
			}
		}
		
		// METHOD 2. check neighbors for closer DP (use DC connectivity)
		/*
		int edge = SphP[sphInd].first_connection;
		int last_edge = SphP[sphInd].last_connection;
		int neighbor = -1;
		
		while(edge >= 0) {
		  neighbor = DC[edge].index;
			
			celldist = Vector(pos.x - P[neighbor].Pos[0],pos.y - P[neighbor].Pos[1]	pos.z - P[neighbor].Pos[2]);
                                           
			if(celldist.PeriodicLengthSquared() < dist2point - INSIDE_EPS)
				sphInd_min[1] = DC[edge].index;
				
			// move to next neighbor
			if(edge == last_edge)
				break;
				
			if (DC[edge].next == edge || DC[edge].next < 0 || dp_neighbor == -5)
			  terminate(" what is going on (%d %d %d) ",DC[edge].next,edge,dp_neighbor);
				
			edge = DC[edge].next;
		}
		*/
		
#ifdef DEBUG
		cout << "VerifyPointInCell PASSED! pt.x = " << pos.x << " pt.y = " << pos.y << " pt.z = " << pos.z << endl;
		
		cout << " [ ]   sphInd_cur = " << parInd << " P.x = " << P[parInd].Pos[0] << " P.y = " << P[parInd].Pos[1] 
				 << " P.z = " << P[parInd].Pos[2] << " (dist2point = " << dist2point << ")" << endl;
#endif // DEBUG

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
	 
#ifdef DEBUG_VERIFY_ENTRY_CELLS
		int next_tetra = -1;
		IF_DEBUG(cout << " TETRA tt = " << tt	<< " moves = " << moves << " ret = " << ret << endl);
		
		ret = -1;
		int check = InTetra(T, tt, &p, &ret, &next_tetra);
		
		if (check < 1) // 2 or 3 ok (on face or edge)
			terminate("ERROR: Entry tetra check failed.");
#endif

		// sanity check
		if (tt < 0 || tt >= Ndt)
			terminate("ERROR: Entry tetra search ended at tt = %d",tt);
		
		ray.tetra = tt;
		*prevEntryTetra = tt;
}

void addValsContribution( vector<float> &vals, int SphP_ind, double weight )
{
		if( weight < INSIDE_EPS ) // catches negative weights as well
			return;
		
		vals[TF_VAL_DENS]    += SphP[SphP_ind].Density * weight;
		vals[TF_VAL_TEMP]    += SphP[SphP_ind].Utherm * weight;
		vals[TF_VAL_VMAG]    += sqrt( P[SphP_ind].Vel[0] * P[SphP_ind].Vel[0] + 
														    	P[SphP_ind].Vel[1] * P[SphP_ind].Vel[1] + 
																	P[SphP_ind].Vel[2] * P[SphP_ind].Vel[2] ) * weight;
		vals[TF_VAL_ENTROPY] += SphP[SphP_ind].OldMass * weight;
		vals[TF_VAL_METAL]   += SphP[SphP_ind].Metallicity * weight;
		// SZ y-parameter (no constants, not real units)
		// line integral of: n_e*T = (x_e*rho*T) = (x_e*rho*(U/mu))
		// require/assume: Utherm is temperature in units of Kelvin already
		vals[TF_VAL_SZY]     += (SphP[SphP_ind].Density * P[SphP_ind].OldAcc) * 
												    SphP[SphP_ind].Utherm * weight;
		// Xray (no constants, not real units): line integral of n^2 * Lambda(T)
		// no additional weighting as with (Temp,Vmag,Ent,Metal)
		// assume: Utherm is temperature in units of Kelvin already (restriction: hot gas only, >1e6)
		vals[TF_VAL_XRAY]    +=  sqrt( SphP[SphP_ind].Utherm ) *
														 pow( SphP[SphP_ind].Density * P[SphP_ind].OldAcc, 2 ) * weight;
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
	cout << "FindNearestGasParticle(): found nearest = " << nearest << " x = " << P[nearest].Pos[0] << " y = " << 
			P[nearest].Pos[1] << " z = " << P[nearest].Pos[2] << " cur_mindist = " << cur_mindist << endl;
#endif	
	
	*mindist = sqrt(cur_mindist_sq);
	
	return nearest;
}

// get primary hydro IND - handle local ghosts
inline int ArepoMesh::getSphPID(int dpInd)
{
		int sphInd = -1;
		
		if (dpInd >= 0 && dpInd < NumGas)
				sphInd = dpInd;
		else if (dpInd >= NumGas)
				sphInd = dpInd - NumGas;
		
		if (sphInd < 0)
				terminate("Negative sphInd in getSphPID().");
				
		return sphInd;
}

void ArepoMesh::checkCurCellTF(bool *addFlag, int sphInd, vector<float> &vals)
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
					vals_next[TF_VAL_DENS] = (float) SphP[qmin].Density;
					vals_next[TF_VAL_TEMP] = (float) SphP[qmin].Utherm;
					
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
	if (fabs(pt.x - P[ray.index].Pos[0]) <= INSIDE_EPS && 
	    fabs(pt.y - P[ray.index].Pos[1]) <= INSIDE_EPS &&
			fabs(pt.z - P[ray.index].Pos[2]) <= INSIDE_EPS)
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
	
#ifdef DEBUG_VERIFY_INCELL_EACH_STEP
	if ( next_tetra != ray.tetra )
	{
		IF_DEBUG(cout << "  TETRA ADVANCE old = " << ray.tetra << " new = " << next_tetra << endl);
				
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
#endif

	ray.tetra = next_tetra;
}

bool ArepoMesh::AdvanceRayOneCellNew(const Ray &ray, double *t0, double *t1, 
																		 Spectrum &Lv, Spectrum &Tr, int threadNum)
{
		double min_t = MAX_REAL_NUMBER;
		int qmin = -1, qmin_dp = -1; // next primary cell SphP/DP index
		
		// verify task
		if (ray.task != ThisTask)
			terminate("Ray on wrong task.");
	
	  Point pos = ray(ray.min_t);

		int SphP_ID = ray.index;
		
		double length; // find_next_voronoi_cell() return
		
		double dir[3];
		dir[0] = ray.d[0];
		dir[1] = ray.d[1];
		dir[2] = ray.d[2];
		
		// TODO: change ray.index to ray.prev_index
		qmin = find_next_cell_DC(T, SphP_ID, &(pos[0]), dir, ray.index, &length);
		
		if( qmin != -1 )
			qmin_dp = DC[qmin].dp_index; // DP_index (in ray.index we store SphP_index)
		
		qmin = DC[qmin].index; // SphP_index (already NumGas subtracted if necessary)
		min_t = ray.min_t + length;
		
#ifdef DEBUG
		cout << "  NEW intersection t = " << min_t << " setting new min_t, next_sphID = " << qmin 
		     << " next_dpID = " << qmin_dp << endl;
#endif

#ifdef DEBUG_VERIFY_INCELL_EACH_STEP
		// verify ray is where we expect it (halfway through this cell)
		Point checkPos = ray( 0.5*(ray.min_t + min_t) );
		ArepoMesh::VerifyPointInCell(ray.index,checkPos);
#endif
	
		// check if exiting box and failed to exit a face
		if (qmin == -1) {
			Point exitcell = ray(min_t);
			
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
		if (qmin != -1 && SphP_ID >= NumGas)
			addFlag = false;
		
		// pack cell-center values to test TF
		vector<float> vals( TF_NUM_VALS, 0.0 );
		
		addValsContribution( vals, SphP_ID, 1.0 );
		
		// if TF evaluates to zero at this cell midpoint, and at all neighboring cell midpoints, 
		// then we are guaranteed that nowhere in this cell will be important, so jump onwards
		checkCurCellTF(&addFlag, SphP_ID, vals);

		// check for proper exit point
		if (qmin != -1)
		{
			// clamp min_t to avoid integrating outside the box
			IF_DEBUG(cout << " have exit: min_t = " << min_t << " (t1=" << *t1 << " t0=" << *t0 << ") addFlag = " << addFlag << endl);
			min_t = Clamp(min_t,*t0,*t1);
			
			if (addFlag)
			{
				// entry and exit points for this cell
				Point hitcell  = ray(ray.min_t);
				Point exitcell = ray(min_t);
					
				// cell gradient information
				Vector sphCen(     SphP[SphP_ID].Center[0],   SphP[SphP_ID].Center[1],   SphP[SphP_ID].Center[2]);
				Vector sphDensGrad(SphP[SphP_ID].Grad.drho[0],SphP[SphP_ID].Grad.drho[1],SphP[SphP_ID].Grad.drho[2]);
									
				IF_DEBUG(hitcell.print("  hcell "));
				IF_DEBUG(exitcell.print("  ecell "));
				IF_DEBUG(sphCen.print("  dpCen "));
									
				Vector norm = exitcell - hitcell;
					
				// compute total path length through cell
				double len = norm.Length();
					
				// find interpolated density at midpoint of line segment through voronoi cell
				Vector midpt(hitcell[0] + 0.5 * norm[0],hitcell[1] + 0.5 * norm[1],hitcell[2] + 0.5 * norm[2]);
				midpt -= sphCen;

				// optical depth: always use gradient for tau calculation (though apply as one value for whole cell)
				Spectrum stepTau(0.0);
				stepTau += transferFunction->sigma_t() * (SphP[SphP_ID].Density + Dot(sphDensGrad,midpt)) * len;
				//Tr *= Exp(-stepTau); // reduce transmittance for optical depth
					
				// setup sampling for this cell
				Point prev_sample_pt;
				int nSamples;
				float stepSize;
				
				if( Config.viStepSize > 0.0 )
				{
					stepSize = Config.viStepSize;
					
					// sub-stepping: strict in world space
					prev_sample_pt = ray(*t0 + ray.depth * stepSize);
					norm = exitcell - prev_sample_pt;

					nSamples = (int)floor(norm.Length() / stepSize);
					norm = Normalize(norm);
				}
				else if( Config.viStepSize < 0.0 )
				{
					// sub-stepping: adaptive based on cell size (could include gradients as well)
					// in this case -Config.viStepSize is the NUMBER of samples per cell
					// note: Config.viStepSize=-1 should reproduce the behavior of Config.viStepSize=0
					stepSize = len / -Config.viStepSize;
					nSamples = (int)(-Config.viStepSize); // = floor(len/stepSize) up to floating rounding
					norm = Normalize(norm);
					prev_sample_pt = hitcell - 0.5 * stepSize * norm; // all samples interior to bounding faces
				}
				else
				{
					// not sub-stepping, then do single sample at midpoint of line through cell
					nSamples = 1;
					stepSize = len;
					norm = Normalize(norm);
					prev_sample_pt = hitcell - 0.5 * stepSize * norm;
				}
				
				IF_DEBUG(prev_sample_pt.print("  prev_sample_pt "));
					
				IF_DEBUG(cout << " sub-stepping len = " << len << " nSamples = " << nSamples 
											<< " (step = " << len/nSamples << ")" << endl);
												
				for (int i = 0; i < nSamples; ++i)
				{
					// where are we inside the cell?
					Point samplept = prev_sample_pt + (i+1)*stepSize * norm;
					
#if defined(DTFE_INTERP) || defined(NNI_WATSON_SAMBRIDGE) || defined(NNI_LIANG_HALE)
					locateCurrentTetra(ray, samplept);
#endif

					// temporary fix for near-box-edge rays going crazy with NNI
					if (!extent.Inside(samplept))
					  terminate("ERROR: Sample point outside box. (%g %g %g)",samplept.x,samplept.y,samplept.z);
																							
					// subsample (replace fields in vals by interpolated values)
					int status = subSampleCell(ray, samplept, vals, threadNum);
							
#ifdef DEBUG
					double fracstep = 1.0 / nSamples;
					cout << "  segment[" << i << "] fractrange [" << (i*fracstep) << "," 
							<< (i*fracstep)+fracstep << "] rho = " << SphP[SphP_ID].Density
							<< " rho subSample = " << vals[TF_VAL_DENS];
					samplept.print(" ");
#endif
					
					// compute emission-only source term using transfer function
					if(status)
					  Lv += Tr * transferFunction->Lve(vals) * stepSize;
							
					// update raw column density integrals
					if( Config.projColDens ) {
						// ensure positivity of integral weights
						if( vals[TF_VAL_DENS] < 0.0 )
							vals[TF_VAL_DENS] = 0.0;
							
						// same weighting procedure as in voronoi_makeimage_new()
						// e.g. weight (Temp,Vmag,Ent,Metal) by rho*len and then normalize out at the end
						float weight = vals[TF_VAL_DENS] * stepSize;
						
						// column density
						ray.raw_vals[0] += vals[TF_VAL_DENS] * stepSize;
						
						// mass-weighted values
						ray.raw_vals[1] += vals[TF_VAL_TEMP] * weight;
						ray.raw_vals[2] += vals[TF_VAL_VMAG] * weight;
						ray.raw_vals[3] += vals[TF_VAL_ENTROPY] * weight;
						ray.raw_vals[4] += vals[TF_VAL_METAL] * weight;
						
						// un-weighted line integrals
						ray.raw_vals[5] += vals[TF_VAL_SZY] * stepSize;
						if( vals[TF_VAL_TEMP] >= 1e6 ) // xray restriction: hot gas only (Temp > 1e6 K)
							ray.raw_vals[6] += vals[TF_VAL_XRAY] * stepSize;
					}
							
					// update previous sample point marker
					ray.depth++;
				} // nSamples
			} //addFlag
			
			// update ray: transfer to next voronoi cell (possibly on different task)
			ray.task  = DP[qmin_dp].task;
			ray.prev_index = ray.index;
			ray.index = qmin;
			ray.min_t = Clamp(min_t,ray.min_t,ray.max_t);
			
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
					cout << " P[ray.index] pos: " << P[ray.index].Pos[0] << " " << P[ray.index].Pos[1]
							 << " " << P[ray.index].Pos[2] << endl << endl;
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

void ArepoMesh::ComputeQuantityBounds()
{
		float pmax  = -INFINITY;
		float pmin  = INFINITY;
		float pmean = 0.0;
		
		float umax = -INFINITY;
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
		
		valBounds[TF_VAL_TEMP*3 + 0] = umin;
		valBounds[TF_VAL_TEMP*3 + 1] = umax;
		valBounds[TF_VAL_TEMP*3 + 2] = umean / NumGas;
		
		cout << " Density min = " << valBounds[TF_VAL_DENS*3 + 0] 
									<< " max = " << valBounds[TF_VAL_DENS*3 + 1] 
									<< " mean = " << valBounds[TF_VAL_DENS*3 + 2] << endl;
									
		cout << " Temp  min = " << valBounds[TF_VAL_TEMP*3 + 0] 
							 << " max = " << valBounds[TF_VAL_TEMP*3 + 1] 
							 << " mean = " << valBounds[TF_VAL_TEMP*3 + 2] << endl;
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
		//	cout << "[" << i << "] numVert=" << numVertices[i-1] << " offset=" << vertexOffset.back() << endl;
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
														 << " oldmass/entr = " << SphP[i].OldMass << endl;		
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
				cout << " SphP[" << setw(3) << i << "] DC.first = " << setw(2) << SphP[i].first_connection
				     << " (" << DC[SphP[i].first_connection].index << ")";
				
				do {
						next = DC[dc].next;
						cout << " " << " next = " << setw(2) << next << " (" << DC[next].index << ")";
						dc = next;
				} while (next != SphP[i].last_connection);
				
				cout << " DC.last = " << SphP[i].last_connection 
				     << " (" << DC[SphP[i].last_connection].index << ")" << endl;		
		}
		
		if (numVertices.size()) {
			cout << endl << "Voronoi Edges (NumGas=" << NumGas << "):" << endl;
			for (int i=0; i < NumGas; i++) {
					cout << setw(3) << i << " numVert = " << setw(2) << numVertices[i] << " vertexOffset = " 
							 << setw(3) << vertexOffset[i] << endl;
			}
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
