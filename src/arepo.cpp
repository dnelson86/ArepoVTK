/*
 * arepo.cpp
 * dnelson
 */
 
#include <alloca.h>
 
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
		
		// TODO: temp units
		unitConversions[TF_VAL_DENS]   = All.UnitDensity_in_cgs / MSUN_PER_PC3_IN_CGS;
		unitConversions[TF_VAL_UTHERM] = All.UnitEnergy_in_cgs;
		
		IF_DEBUG(cout << "unitConv[dens]   = " << unitConversions[TF_VAL_DENS] << endl);
		IF_DEBUG(cout << "unitConv[utherm] = " << unitConversions[TF_VAL_UTHERM] << endl);
		
		// debugging
		//ArepoMesh::DumpMesh();
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
}

void ArepoMesh::setupAuxMeshes(void)
{
#ifdef NATURAL_NEIGHBOR_INTERP
		// allocate auxiliary meshes
		int numMeshes = Config.nTasks;
		AuxMeshes = new tessellation[numMeshes];
		
		// define neutral index
		DPinfinity = -5;		
		
		double box = All.BoxSize;
		double tetra_incircle, tetra_sidelength, tetra_height, tetra_face_height;		
		
		tetra_incircle = 1.5 * box;
		tetra_sidelength = tetra_incircle * sqrt(24);
		tetra_height = sqrt(2.0 / 3) * tetra_sidelength;
		tetra_face_height = sqrt(3.0) / 2.0 * tetra_sidelength;		
		
		for(int k = 0; k < numMeshes; k++)
		{
			// copy individual allocation factors for auxiliary meshes
			//init_clear_auxmesh(&AuxMeshes[k]);

			point *p;
			int i;
			
			AuxMeshes[k].Indi.AllocFacNdp = AUXMESH_ALLOC_SIZE / 2;
			AuxMeshes[k].Indi.AllocFacNdt = AUXMESH_ALLOC_SIZE;
			AuxMeshes[k].Indi.AllocFacNvf = AUXMESH_ALLOC_SIZE;

			AuxMeshes[k].MaxNdp = (int)T->Indi.AllocFacNdp;
			AuxMeshes[k].MaxNdt = (int)T->Indi.AllocFacNdt;
			AuxMeshes[k].MaxNvf = (int)T->Indi.AllocFacNvf;		

			AuxMeshes[k].VF = static_cast<face*>
										 (mymalloc_movable(AuxMeshes[k].VF, "VFaux", AuxMeshes[k].MaxNvf * sizeof(face)));
			AuxMeshes[k].DP = static_cast<point*>
										 (mymalloc_movable(AuxMeshes[k].DP, "DPaux", (AuxMeshes[k].MaxNdp+5) * sizeof(point)));
			AuxMeshes[k].DP += 5; // leave first five for bounding tetra + infinity
			AuxMeshes[k].DT = static_cast<tetra*>
										 (mymalloc_movable(AuxMeshes[k].DT, "DTaux", AuxMeshes[k].MaxNdt * sizeof(tetra)));
			
			AuxMeshes[k].DTC = static_cast<tetra_center*>(
											mymalloc_movable(&AuxMeshes[k].DTC, "AuxDTC", AuxMeshes[k].MaxNdt * sizeof(tetra_center)));
			AuxMeshes[k].DTF = static_cast<char*>(
											mymalloc_movable(&AuxMeshes[k].DTF, "AuxDTF", AuxMeshes[k].MaxNdt * sizeof(char)));	

			point *DP = AuxMeshes[k].DP;
			//tetra *DT = AuxMeshes[k].DT;

			/* first, let's make the points */
			DP[-4].x = 0.5 * tetra_sidelength;
			DP[-4].y = -1.0 / 3 * tetra_face_height;
			DP[-4].z = -0.25 * tetra_height;

			DP[-3].x = 0;
			DP[-3].y = 2.0 / 3 * tetra_face_height;
			DP[-3].z = -0.25 * tetra_height;

			DP[-2].x = -0.5 * tetra_sidelength;
			DP[-2].y = -1.0 / 3 * tetra_face_height;
			DP[-2].z = -0.25 * tetra_height;

			DP[-1].x = 0;
			DP[-1].y = 0;
			DP[-1].z = 0.75 * tetra_height;

			for(i = -4; i <= -1; i++)
				{
					DP[i].x += 0.5 * box;
					DP[i].y += 0.5 * box;
					DP[i].z += 0.5 * box;
				}

			for(i = -4, p = &DP[-4]; i < 0; i++, p++)
				{
					p->index = -1;
					p->task = ThisTask;
					p->timebin = 0;
				}

			/* we also define a neutral element at infinity */
			DP[DPinfinity].x = GSL_POSINF;
			DP[DPinfinity].y = GSL_POSINF;
			DP[DPinfinity].z = GSL_POSINF;
			DP[DPinfinity].index = -1;
			DP[DPinfinity].task = ThisTask;
			DP[DPinfinity].timebin = 0;
			
		}
#endif
}

void ArepoMesh::precomputeTetraGrads()
{
#ifdef DTFE_INTERP
		int i;
		double det_A;
		double A[3][3]; // row,col
		double A_inv[3][3];
		float delta_f[3];
		float field_0;
		
		// allocate for grad[3] per tetra (one field only)
		DT_grad = new float[3 * Ndt];
		
		for (i = 0; i < Ndt; i++) {
				// if this tetra is connected to the bounding tetra then skip it
				if ( DT[i].p[0] < 0 || DT[i].p[1] < 0 || DT[i].p[2] < 0 || DT[i].p[3] < 0 )
				  continue;
					
				// fill A
				A[0][0] = DP[DT[i].p[1]].x - DP[DT[i].p[0]].x;
				A[1][0] = DP[DT[i].p[2]].x - DP[DT[i].p[0]].x;
				A[2][0] = DP[DT[i].p[3]].x - DP[DT[i].p[0]].x;
				A[0][1] = DP[DT[i].p[1]].y - DP[DT[i].p[0]].y;
				A[1][1] = DP[DT[i].p[2]].y - DP[DT[i].p[0]].y;
				A[2][1] = DP[DT[i].p[3]].y - DP[DT[i].p[0]].y;
				A[0][2] = DP[DT[i].p[1]].z - DP[DT[i].p[0]].z;
				A[1][2] = DP[DT[i].p[2]].z - DP[DT[i].p[0]].z;
				A[2][2] = DP[DT[i].p[3]].z - DP[DT[i].p[0]].z;
				
				// compute determinant(A)
				det_A = A[0][0]*A[1][1]*A[2][2] + A[1][0]*A[2][1]*A[0][2] +
                A[2][0]*A[0][1]*A[1][2] - A[0][0]*A[2][1]*A[1][2] - 
								A[2][0]*A[1][1]*A[0][2] - A[1][0]*A[0][1]*A[2][2];
								
				if (det_A < INSIDE_EPS && det_A > -INSIDE_EPS) { // debugging
				  cout << "ERROR: det_A is zero! " << det_A << " [" << i << "]" << endl;
					exit(12509);
				}
				
				// reciprocal of determinant
				det_A = 1.0 / det_A;
				
				// compute A_inv
				A_inv[0][0] = det_A * ( A[1][1]*A[2][2] - A[1][2]*A[2][1] );
				A_inv[1][0] = det_A * ( A[1][2]*A[2][0] - A[1][0]*A[2][2] );
				A_inv[2][0] = det_A * ( A[1][0]*A[2][1] - A[1][1]*A[2][0] );
				A_inv[0][1] = det_A * ( A[0][2]*A[2][1] - A[0][1]*A[2][2] );
				A_inv[1][1] = det_A * ( A[0][0]*A[2][2] - A[0][2]*A[2][0] );
				A_inv[2][1] = det_A * ( A[0][1]*A[2][0] - A[0][0]*A[2][1] );
				A_inv[0][2] = det_A * ( A[0][1]*A[1][2] - A[0][2]*A[1][1] );
				A_inv[1][2] = det_A * ( A[0][2]*A[1][0] - A[0][0]*A[1][2] );
				A_inv[2][2] = det_A * ( A[0][0]*A[1][1] - A[0][1]*A[1][0] );
				
				// compute delta_f (field specific)
				field_0 = SphP[ getSphPID(DP[DT[i].p[0]].index) ].Density;
				
				delta_f[0] = SphP[ getSphPID(DP[DT[i].p[1]].index) ].Density - field_0;
				delta_f[1] = SphP[ getSphPID(DP[DT[i].p[2]].index) ].Density - field_0;
				delta_f[2] = SphP[ getSphPID(DP[DT[i].p[3]].index) ].Density - field_0;
				
				// compute the gradient (A_inv * delta_f) and store x,y,z components in DT_grad
				DT_grad[3*i+0] = A_inv[0][0]*delta_f[0] + A_inv[0][1]*delta_f[1] + A_inv[0][2]*delta_f[2];
				DT_grad[3*i+1] = A_inv[1][0]*delta_f[0] + A_inv[1][1]*delta_f[1] + A_inv[1][2]*delta_f[2];
				DT_grad[3*i+2] = A_inv[2][0]*delta_f[0] + A_inv[2][1]*delta_f[1] + A_inv[2][2]*delta_f[2];
		}
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
		
		// exactly on DP point, avoid divide by zero
		if (dist2point == 0.0)
				return;
				
		double dist2point_min = dist2point;
		
		// check neighbors for a closer tetra midpoint
		const int start_edge = midpoint_idx[dp].first;
		const int n_edges    = midpoint_idx[dp].second;
		
		for (int i=0; i < n_edges; i++)
		{
				const int dp_neighbor = opposite_points[start_edge + i];
				
				celldist = Vector(pos.x - DP[dp_neighbor].x,
				                  pos.y - DP[dp_neighbor].y,
													pos.z - DP[dp_neighbor].z);
				
				dist2point = celldist.LengthSquared();
				
				if (dist2point < dist2point_min) {
						dist2point_min = dist2point;
				}
		}
	
		// check
		if (dist2point/dist2point_min > 1+INSIDE_EPS) {
				cout << "VerifyPointInCell FAILED! dp = " << dp << " pos.x = " << pos.x << " pos.y = " << pos.y << " pos.z = " << pos.z << endl;
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
	//if(image_get_next_tetra(T, ray.tetra, &pp, &pend, &next_tetra, &pnew, &previous_tetra))
	if ( next_tetra != ray.tetra )
	{
#ifdef DEBUG
		cout << "  TETRA ADVANCE old = " << ray.tetra << " new = " << next_tetra << endl;
				
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
		
		if(test > 1)
			cout << "  TETRA: WARNING: on face or edge [" << test << "]" << endl;
#endif
		//previous_tetra = ray.tetra;
		ray.tetra = next_tetra;
	}
			
	// update pp (pstart of the line) to this midpt
	//pp = pend;
}

bool ArepoMesh::AdvanceRayOneCellNew(const Ray &ray, float *t0, float *t1, 
																		 int previous_cell, Spectrum &Lv, Spectrum &Tr, int taskNum)
{
		double min_t = MAX_REAL_NUMBER;
		int qmin = -1; // next
		
		// verify task
		if (ray.task != ThisTask) {
				cout << "[" << ThisTask << "] ERROR! ray.task = " << ray.task << " differs from ThisTask!" << endl;
				terminate("1138");
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
		
#ifdef DEBUG
		Point exitbox = ray(*t1);
		
		hitbox.print(" hitbox ");
		exitbox.print(" exitbox ");
		cellp.print(" cellp (DP.xyz) ");
#endif

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
	
			IF_DEBUG(midp.print("  midp (midpoints-hitbox) "));
			IF_DEBUG(norm.print("  norm (midpoints-cellp) "));
	
			// find intersection of ray with this face
			double dotprod1 = Dot( ray.d, norm ); //ddotq
			double dotprod2 = Dot( midp,  norm ); //cdotq
			
			// check if ray is aligned on face (e.g. backgroundgrid)
			if (dotprod1 == 0 && dotprod2 == 0)
				continue;
			
			if (dotprod1 > 0)
			{
				double t = dotprod2 / dotprod1; // infinite line/plane intersection test
				
				IF_DEBUG(cout << "  i[" << i << "] dotprod>0 and t = " << t << " (min=" << (ray.min_t-*t0) << ")" << endl);
				
				if (t > (ray.min_t-*t0) && t < min_t) {
					min_t = t;
					qmin = opposite_points[edge.first + i];
					IF_DEBUG(cout << "  intersection t = " << t << " setting new min_t, qmin (next DP) = " << qmin << endl);
				}
		} else {
				//cout << "Warning: Point on wrong side of face." << endl;
				// check if numerical error has put it on the wrong side of the face
				if (dotprod2 > 0) {
					//terminate("1137");
					//TODO: check and allow if we are really close to a face
				}
			}
		}
		
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
			
#if defined(DTFE_INTERP) || defined(NATURAL_NEIGHBOR_WATSON)
				//int previous_tetra = ray.tetra;

				//pp.x = hitcell.x;
				//pp.y = hitcell.y;
				//pp.z = hitcell.z;
				//set_integers_for_pointer(&pp);
					
				//pend.x = midpt.x;
				//pend.y = midpt.y;
				//pend.z = midpt.z;
							
				//set_integers_for_pointer(&pend);
					
				//point pnew; // image_get_next_tetra return
#endif

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
					
#if defined(DTFE_INTERP) || defined(NATURAL_NEIGHBOR_WATSON)
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
					terminate("1130");
			}
		} // qmin
		
		return true;		
}

#ifdef NATURAL_NEIGHBOR_SPHKERNEL
inline float sph_kernel(float dist, float hinv)
{
  float u = dist * hinv;
  
  if(u < 0.5)
    return (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
  if(u < 1.0)
    return KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
  return 0;
}
#endif

// interpolate scalar fields at position pt inside Voronoi cell SphP_ID (various methods)
int ArepoMesh::subSampleCell(int SphP_ID, const Ray &ray, Vector &pt, float *vals, int taskNum)
{
	int DP_ID = ray.index;
	
	// check degenerate point in R3, immediate return
	if (fabs(pt.x - DP[DP_ID].x) <= INSIDE_EPS &&
			fabs(pt.y - DP[DP_ID].y) <= INSIDE_EPS &&
			fabs(pt.z - DP[DP_ID].z) <= INSIDE_EPS)
	{
			vals[TF_VAL_DENS]   = SphP[SphP_ID].Density;
			vals[TF_VAL_UTHERM] = SphP[SphP_ID].Utherm;
			return 1;
	}
				
	// zero vals we will override in this function
	vals[TF_VAL_DENS]   = 0.0;
	vals[TF_VAL_UTHERM] = 0.0;			

  // for periodic distances
	float dx,dy,dz,xtmp,ytmp,ztmp;	
				
#ifdef NATURAL_NEIGHBOR_IDW

/* exponent of distance, greater values assign more influence to values closest to the 
 * interpolating point, approaching piecewise constant for large POWER_PARAM.
 * in N dimensions, if p <= N, the interpolated values are dominated by points far away,
 * which is rather bizarre.
 */
 
#define POWER_PARAM 4.0

		// list of neighbors
		const int start_edge = midpoint_idx[DP_ID].first;
		const int n_edges    = midpoint_idx[DP_ID].second;
		
		// add parent to list
		int dp_neighbor, sphp_neighbor;
		float weight,weightsum=0,distsq;
		
		// loop over each neighbor
		for (int k=0; k < n_edges; k++) {
			dp_neighbor = opposite_points[start_edge + k];
			sphp_neighbor = getSphPID(dp_neighbor);
			
			// calculate weight as 1.0/dist^power (Shepard's Method)
			dx = NGB_PERIODIC_LONG_X(DP[dp_neighbor].x - pt.x);
			dy = NGB_PERIODIC_LONG_Y(DP[dp_neighbor].y - pt.y);
			dz = NGB_PERIODIC_LONG_Z(DP[dp_neighbor].z - pt.z);
			
			distsq = dx*dx + dy*dy + dz*dz;
			
			weight = 1.0 / pow(distsq,POWER_PARAM);
			weightsum += weight;
			
			vals[TF_VAL_DENS]   += SphP[sphp_neighbor].Density * weight;
			vals[TF_VAL_UTHERM] += SphP[sphp_neighbor].Utherm * weight;
		}
		
		// add in primary parent
		dx = NGB_PERIODIC_LONG_X(DP[DP_ID].x - pt.x);
		dy = NGB_PERIODIC_LONG_Y(DP[DP_ID].y - pt.y);
		dz = NGB_PERIODIC_LONG_Z(DP[DP_ID].z - pt.z);
	
		distsq = dx*dx + dy*dy + dz*dz;
						 
		weight = 1.0 / pow(distsq,POWER_PARAM);
		weightsum += weight;
		
		vals[TF_VAL_DENS]   += SphP[SphP_ID].Density * weight;
		vals[TF_VAL_UTHERM] += SphP[SphP_ID].Utherm * weight;
		
		// normalize weights
		vals[TF_VAL_DENS]   /= weightsum;
		vals[TF_VAL_UTHERM] /= weightsum;
			
#endif // NATURAL_NEIGHBOR_IDW

#ifdef NATURAL_NEIGHBOR_SPHKERNEL

/* use <1 for more smoothing, makes hsml_used bigger than hsml_ngb_max
 * use >1 for less smoothing, make hsml_used smaller than hsml_ngb_max 
 *   (at least some some neighbors will not contribute, and if this is too big, the
 *    containing cell may also not contribute, leading to black holes)
 */
#define HSML_FAC 0.95

		// list of neighbors
		const int start_edge = midpoint_idx[DP_ID].first;
		const int n_edges    = midpoint_idx[DP_ID].second;
		
		// add parent to list
		int dp_neighbor, sphp_neighbor;
		float weight,weightsum=0,distsq;
		
    // 1. smoothing length parameter: calculate over neighbor distances
    float hsml2 = 0.0;
    for (int k=0; k < n_edges; k++) {
		  dp_neighbor = opposite_points[start_edge + k];
			
	    if(dp_neighbor < 0)
        continue;

			dx = NGB_PERIODIC_LONG_X(DP[dp_neighbor].x - pt.x);
			dy = NGB_PERIODIC_LONG_Y(DP[dp_neighbor].y - pt.y);
			dz = NGB_PERIODIC_LONG_Z(DP[dp_neighbor].z - pt.z);
			
			distsq = dx*dx + dy*dy + dz*dz;
                                           
      if(distsq > hsml2)
        hsml2 = distsq;
    }
           
    //if(hsml2 < 0.1)
	  //	hsml2 = 0.1;
 
    float hinv = HSML_FAC / sqrtf(hsml2);
            
		// 2. loop over each neighbor and add contribution
		for (int k=0; k < n_edges; k++) {
			dp_neighbor = opposite_points[start_edge + k];
			sphp_neighbor = getSphPID(dp_neighbor);
			
			// calculate weight as sphkernel(distsq/h) (cubic spline kernel)
			dx = NGB_PERIODIC_LONG_X(DP[dp_neighbor].x - pt.x);
			dy = NGB_PERIODIC_LONG_Y(DP[dp_neighbor].y - pt.y);
			dz = NGB_PERIODIC_LONG_Z(DP[dp_neighbor].z - pt.z);
			
			distsq = dx*dx + dy*dy + dz*dz;
			
			weight = sph_kernel(sqrtf(distsq),hinv);
			weightsum += weight; // * SphP[SphP_ID].Volume;
			
			vals[TF_VAL_DENS]   += SphP[sphp_neighbor].Density * weight;
			vals[TF_VAL_UTHERM] += SphP[sphp_neighbor].Utherm * weight;
		}
		
		// add in primary parent
		dx = NGB_PERIODIC_LONG_X(DP[DP_ID].x - pt.x);
		dy = NGB_PERIODIC_LONG_Y(DP[DP_ID].y - pt.y);
		dz = NGB_PERIODIC_LONG_Z(DP[DP_ID].z - pt.z);
		
		distsq = dx*dx + dy*dy + dz*dz;
		
		weight = sph_kernel(sqrtf(distsq),hinv);
		weightsum += weight; // * SphP[SphP_ID].Volume;
		
		vals[TF_VAL_DENS]   += SphP[SphP_ID].Density * weight;
		vals[TF_VAL_UTHERM] += SphP[SphP_ID].Utherm * weight;
		
		// 3. normalize by weight totals
		vals[TF_VAL_DENS]   /= weightsum;
		vals[TF_VAL_UTHERM] /= weightsum;

#endif // NATURAL_NEIGHBOR_SPHKERNEL

#ifdef NATURAL_NEIGHBOR_INTERP
		int tlast = 0;
		int dp_neighbor, sphp_neighbor;

		float weight;	
		
		double dp_old_vol[AUXMESH_ALLOC_SIZE/2];
		double dp_new_vol[AUXMESH_ALLOC_SIZE/2];
		
		// list of neighbors
		const int start_edge = midpoint_idx[DP_ID].first;
		const int n_edges    = midpoint_idx[DP_ID].second;

		// recreate the voronoi cell of the parent's neighbors (auxiliary mesh approach):
		init_clear_auxmesh(&AuxMeshes[taskNum]);
			
		// construct new auxiliary mesh around pt
		for(int k = 0; k < n_edges; k++)
		{
				int q = opposite_points[start_edge + k];
		
				if(AuxMeshes[taskNum].Ndp + 2 >= AuxMeshes[taskNum].MaxNdp)
					terminate("1157");
				
				// insert point
				AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp] = Mesh.DP[q];
				
		    set_integers_for_point(&AuxMeshes[taskNum], AuxMeshes[taskNum].Ndp);
		    tlast = insert_point_new(&AuxMeshes[taskNum], AuxMeshes[taskNum].Ndp, tlast);
		    AuxMeshes[taskNum].Ndp++;
				
				IF_DEBUG(cout << " inserted neighbor k=" << k << " new tlast=" << tlast << " totnum=" << 
								 AuxMeshes[taskNum].Ndp << endl);
		}
		
		// add primary parent as well
		AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp] = Mesh.DP[DP_ID];
		set_integers_for_point(&AuxMeshes[taskNum], AuxMeshes[taskNum].Ndp);
		tlast = insert_point_new(&AuxMeshes[taskNum], AuxMeshes[taskNum].Ndp, tlast);
		AuxMeshes[taskNum].Ndp++;
		
		// compute old circumcircles and volumes
		compute_circumcircles(&AuxMeshes[taskNum]);
		compute_auxmesh_volumes(&AuxMeshes[taskNum], dp_old_vol);
		
#ifdef DEBUG
		cout << " old volumes:";
		for (int k = 0; k < AuxMeshes[taskNum].Ndp; k++) {
				cout << " [k=" << k << "] " << dp_old_vol[k];
		}
		cout << endl;
#endif		
		
		// add interpolate point
		AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp].x = pt.x;
		AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp].y = pt.y;
		AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp].z = pt.z;
		AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp].ID = -1; //unused
		
		set_integers_for_point(&AuxMeshes[taskNum], AuxMeshes[taskNum].Ndp);

		tlast = insert_point_new(&AuxMeshes[taskNum], AuxMeshes[taskNum].Ndp, tlast);
		AuxMeshes[taskNum].Ndp++;
		
		IF_DEBUG(cout << " inserted interpolate new tlast=" << tlast << " totnum=" << AuxMeshes[taskNum].Ndp << endl);
	
		// compute new circumcircles and volumes
		compute_circumcircles(&AuxMeshes[taskNum]);
		compute_auxmesh_volumes(&AuxMeshes[taskNum], dp_new_vol);

#ifdef DEBUG
		cout << " new volumes:";
		for (int k = 0; k < AuxMeshes[taskNum].Ndp; k++) {
				cout << " [k=" << k << "] " << dp_new_vol[k];
		}
		cout << endl;
#endif

		// calculate scalar value based on neighbor values and area fraction weights
		for (int k=0; k < n_edges; k++) {
			dp_neighbor = opposite_points[start_edge + k];
			sphp_neighbor = getSphPID(dp_neighbor);
			
			weight = dp_old_vol[k] - dp_new_vol[k];
			vals[TF_VAL_DENS]   += SphP[sphp_neighbor].Density * weight;
			vals[TF_VAL_UTHERM] += SphP[sphp_neighbor].Utherm * weight;
		}
		
		// add primary parent
		weight = dp_old_vol[n_edges] - dp_new_vol[n_edges];
		
		vals[TF_VAL_DENS]   += SphP[SphP_ID].Density * weight;
		vals[TF_VAL_UTHERM] += SphP[SphP_ID].Utherm * weight;
		
		// normalize by volume of last added cell (around interp point)
		vals[TF_VAL_DENS]   /= dp_new_vol[AuxMeshes[taskNum].Ndp-1];
		vals[TF_VAL_UTHERM] /= dp_new_vol[AuxMeshes[taskNum].Ndp-1];

#endif // NATURAL_NEIGHBOR_INTERP

#ifdef DTFE_INTERP
		// use DT[tt].p[0] as the x_i sample point
		int tt0_DPID   = DT[ray.tetra].p[0];
		
		// if p[0] is part of the bounding tetra (we are on the edge), cannot do DTFE
		// note: the delaunay triangulation is not periodic, otherwise this would be ok
		if (tt0_DPID < 0) {
			IF_DEBUG(cout << "  dtfe: p[0] = " << tt0_DPID << ", skipping!" << endl);
		  return 0;
		}
		
		int tt0_SphPID = getSphPID(tt0_DPID);
		
		Vector p0cen( DP[tt0_DPID].x, DP[tt0_DPID].y, DP[tt0_DPID].z );
		Vector tetraGrad( DT_grad[3*ray.tetra+0], DT_grad[3*ray.tetra+1], DT_grad[3*ray.tetra+2] );
		
		pt -= p0cen; // make relative to p[0] position (not Voronoi center)

		// apply the (linear) gradient to the sampling point
		vals[TF_VAL_DENS]   = SphP[tt0_SphPID].Density + Dot(tetraGrad,pt);
		
		vals[TF_VAL_UTHERM] = SphP[SphP_ID].Utherm; // TODO
		
#ifdef DEBUG
		cout << "  dtfe: tt = " << ray.tetra << " p0 = " << tt0_DPID << " sph0 = " << tt0_SphPID << endl;
		pt.print("  dtfe: relative point ");
		cout << "  dtfe: dens = " << vals[TF_VAL_DENS] << " sphDens = " << SphP[tt0_SphPID].Density << endl;
#endif
		
#endif

#ifdef CELL_GRADIENTS_ONLY
		// old:
		Vector sphCen(     SphP[SphP_ID].Center[0],   SphP[SphP_ID].Center[1],   SphP[SphP_ID].Center[2]);
		Vector sphDensGrad(SphP[SphP_ID].Grad.drho[0],SphP[SphP_ID].Grad.drho[1],SphP[SphP_ID].Grad.drho[2]);
    pt -= sphCen;	// make relative to cell center	
		
		vals[TF_VAL_DENS]   = SphP[SphP_ID].Density + Dot(sphDensGrad,pt);
		vals[TF_VAL_UTHERM] = SphP[SphP_ID].Utherm; // no utherm gradient unless MATERIALS
		
#endif // CELL_GRADIENTS_ONLY

		return 1;
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
						
						if (DP[dp].index >= 0 && DP[dp].index < NumGas)
								SphP_ID = DP[dp].index;
						else if (dp >= NumGas)
								SphP_ID = DP[dp].index - NumGas;
								
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
				int SphP_ID = -1;
				
				if (DP[i].index >= 0 && DP[i].index < NumGas)
						SphP_ID = DP[i].index;
				else if (i >= NumGas)
						SphP_ID = DP[i].index - NumGas;
				
				if (SphP_ID < 0)
						terminate("1133");
				
				mm.insert(make_pair(SphP_ID, i));
		}
		
		// set up mapping of DP id -> DP primary id
		for (int i=0; i < Ndp; i++)
		{
				int SphP_ID = -1;
				
				if (DP[i].index >= 0 && DP[i].index < NumGas)
						SphP_ID = DP[i].index;
				else if (i >= NumGas)
						SphP_ID = DP[i].index - NumGas;
		
				if (SphP_ID < 0)
						terminate("1134");
		
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
				int SphP_ID = -1;
				
				if (DP[i].index >= 0 && DP[i].index < NumGas)
						SphP_ID = DP[i].index;
				else if (i >= NumGas)
						SphP_ID = DP[i].index - NumGas;
						
				if (SphP_ID < 0)
						continue;
						
				const Vector cellp(DP[i].x,DP[i].y,DP[i].z);
				
				// find connections for this cell
				pair<mmi,mmi> dp_neighbors(conn.equal_range(i));
				
				if (dp_neighbors.first == dp_neighbors.second)
						terminate("1133");
						
				for (; dp_neighbors.first != dp_neighbors.second; dp_neighbors.first++)
				{
						const int dp_neighbor = dp_neighbors.first->second;
				
						int SphP_ID_n = -1;
						
						if (DP[dp_neighbor].index >= 0 && DP[dp_neighbor].index < NumGas)
								SphP_ID_n = DP[dp_neighbor].index;
						else if (dp_neighbor >= NumGas) {
								//cout << " NumGas=" << NumGas << " dp_neighbor=" << dp_neighbor << " DP[dp_neighbor].index="
								//     << DP[dp_neighbor].index << endl;
								SphP_ID_n = DP[dp_neighbor].index - NumGas;
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
