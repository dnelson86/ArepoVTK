/*
 * arepo.cpp
 * dnelson
 */
 
#include "arepo.h"

#ifdef ENABLE_AREPO

// check for required Arepo compilation options

#ifndef VORONOI_NEW_IMAGE
#error ERROR. Missing required Arepo compilation option VORONOI_NEW_IMAGE.
#endif
#ifndef PEANOHILBERT
#error ERROR. Missing required Arepo compilation option PEANOHILBERT.
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
}

bool Arepo::LoadSnapshot()
{
    IF_DEBUG(cout << "Arepo::LoadSnapshot(" << snapFilename << ")." << endl);
		
		freopen("/dev/null","w",stdout); //hide arepo stdout
		
		// set startup options
		WriteMiscFiles = 0;
		RestartSnapNum = -1;
		RestartFlag    = SUNRISE;

		strcpy(ParameterFile,paramFilename.c_str());

		// call arepo: run setup
		begrun1();
		
		// load snapshot
		if (snapFilename.substr(snapFilename.size()-5,5) == ".hdf5")
		  snapFilename = snapFilename.substr(0,snapFilename.size()-5);
			
		read_ic(snapFilename.c_str());
		
		// call arepo: read snapshot, allocate storage for tree, 
		//             initialize particle data, domain decomposition, initial HSML
  	if (init() != SUNRISE) {
				cout << "Arepo::LoadSnapshot() ERROR: Arepo did not return successfully." << endl;
				return false;
		}
	
		// units
		// TODO
	
		freopen("/dev/tty","w",stdout); //return stdout to terminal
		
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
		
		// preprocessing
		//ArepoMesh::ComputeVoronoiEdges();
		ArepoMesh::ComputeQuantityBounds();
		
		// boxsize
		extent = BBox(Point(0.0,0.0,0.0),Point(boxSize,boxSize,boxSize));
		
		IF_DEBUG(extent.print(" ArepoMesh extent "));
		
		// debugging
		//arepoMesh->DumpMesh();
}

void ArepoMesh::LocateEntryCellBrute(const Ray &ray, float *t0, float *t1)
{
		// note: using the brute force search is O(N_rays * NumPart) - not good
		Point hitbox  = ray(*t0);
		
		int corig;
		double mindist = MAX_REAL_NUMBER;
		
		for (int i=0; i < N_gas; i++) {
				double dist = sqrt( (hitbox.x-P[i].Pos[0])*(hitbox.x-P[i].Pos[0]) + 
													  (hitbox.y-P[i].Pos[1])*(hitbox.y-P[i].Pos[1]) + 
														(hitbox.z-P[i].Pos[2])*(hitbox.z-P[i].Pos[2]) );
				if (dist < mindist) {
						mindist = dist;
						corig = i;
				}
		}

		IF_DEBUG(cout << " brute corig = " << corig << " P[corig].index = none" 
									<< " dist = " << mindist << " (x = " << P[corig].Pos[0] << " y = "
									<< P[corig].Pos[1] << " z = " << P[corig].Pos[2] << ")" << endl);
		
		ray.index = corig;
		ray.task  = 0; //TODO
		
		// check for local ghosts
		if (ray.task == ThisTask && ray.index >= N_gas) {
				cout << "ERROR! ray.index = " << ray.index << " (N_gas=" << N_gas << " Ndp=" << Ndp 
						 << ") [task=" << ray.task << "] out of bounds." << endl;
				cout << " hitbox: " << hitbox.x << " " << hitbox.y << " " << hitbox.z << endl;
				cout << " DP nearest pos: " << DP[ray.index].x << " " << DP[ray.index].y << " " << DP[ray.index].z << endl;
				cout << " local N_gas pos: " << P[ray.index-N_gas].Pos[0] << " " << P[ray.index-N_gas].Pos[1]
						 << " " << P[ray.index-N_gas].Pos[2] << endl;
				endrun(1107);
		}
		// verify task assignment
		if (ray.task < 0 || ray.task >= NTask) {
				cout << "ERROR! ray has bad task=" << ray.task << endl;
				endrun(1115);
		}
		
}

void ArepoMesh::LocateEntryCell(const Ray &ray, float *t0, float *t1)
{
		Point hitbox  = ray(*t0);
		Point exitbox = ray(*t1);
		
		IF_DEBUG(cout << " ray hits box at x = " << hitbox.x << " y = " << hitbox.y << " z = " << hitbox.z << endl);
		IF_DEBUG(cout << " ray exit box at x = " << exitbox.x << " y = " << exitbox.y << " z = " << exitbox.z << endl);		
		
		// use peanokey to find domain and task for ray
		if (ray.task == -1) {
				// TODO
				ray.task = 0;
		}
		
		// TODO: exchange
	
		double pos[3], mindist, newdist;
		pos[0] = hitbox.x;
		pos[1] = hitbox.y;
		pos[2] = hitbox.z;
		
		// use tree to find nearest gas particle (local only)
		const int corig = ArepoMesh::FindNearestGasParticle(hitbox, &mindist);
	 
		IF_DEBUG(cout << " corig = " << corig 
									<< " dist = " << mindist << " (x = " << P[corig].Pos[0] << " y = "
									<< P[corig].Pos[1] << " z = " << P[corig].Pos[2] << ")" << endl); 
	 
		// refine nearest point search to account for local ghosts
	/*
		int count = 0;
		int corig_min, dp;
	
		int q = SphP[corig].first_connection;
		
		// look at all neighbors of the starting point (across its voronoi faces)
		while (q >= 0)
		{
				dp = DC[q].dp_index;
				
				// calculate distance
				newdist = sqrt( (pos[0] - P[dp].Pos[0])*(pos[0] - P[dp].Pos[0]) + 
												(pos[1] - P[dp].Pos[1])*(pos[1] - P[dp].Pos[1]) + 
												(pos[2] - P[dp].Pos[2])*(pos[2] - P[dp].Pos[2]) );
				
				cout << " q=" << q << " (neighbor dp=" << dp << ") newdist = " << newdist << endl;
				
				// set new mindist if found
				if (newdist < mindist) {
						cout << " newdist<mindist, setting new corig_min=" << dp << endl;
						mindist = newdist;
						corig_min = dp;
				}
				
				// move to next face
				if (q == SphP[corig].last_connection) {
						cout << " done with faces, breaking" << endl;
						break;
				}
						
				q = DC[q].next;
				count++;
		}	
		
		IF_DEBUG(cout << " iterates count = " << count << " picked ray.index = " << ray.index << endl);

		IF_DEBUG(cout << " corig_min = " << corig_min
									<< " dist = " << mindist << " (x = " << P[corig_min].Pos[0] << " y = "
									<< P[corig_min].Pos[1] << " z = " << P[corig_min].Pos[2] << ")" << endl);
		
		*/
		
		ray.index = corig; //corig_min
		ray.task  = 0; //TODO
							
		// check for local ghosts
		if (ray.task == ThisTask && ray.index >= N_gas) {
				cout << "ERROR! ray.index = " << ray.index << " (N_gas=" << N_gas << " Ndp=" << Ndp 
						 << ") [task=" << ray.task << "] out of bounds." << endl;
				cout << " hit pos: " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
				cout << " DP nearest pos: " << DP[ray.index].x << " " << DP[ray.index].y << " " << DP[ray.index].z << endl;
				cout << " local N_gas pos: " << P[ray.index-N_gas].Pos[0] << " " << P[ray.index-N_gas].Pos[1]
						 << " " << P[ray.index-N_gas].Pos[2] << endl;
				endrun(1107);
		}
		// verify task assignment
		if (ray.task < 0 || ray.task >= NTask) {
				cout << "ERROR! ray has bad task=" << ray.task << endl;
				endrun(1115);
		}
		
}

int ArepoMesh::FindNearestGasParticle(Point &pt, double *mindist)
{
		// based on ngb.c:treefind_nearest_local() but without periodic modifiers
		
		int node, nearest, p;
		struct NODE *current;
		double dx, dy, dz, cur_mindist;

		// top node
		node = All.MaxPart;
		
		// pick random gas particle for guess of the min distance (why?)
		//nearest = floor(get_random_number(SelRnd++) * N_gas);
		nearest = (int)floor(N_gas/2.0);
		dx = P[nearest].Pos[0] - pt.x;
		dy = P[nearest].Pos[1] - pt.y;
		dz = P[nearest].Pos[2] - pt.z;
		cur_mindist = sqrt(dx * dx + dy * dy + dz * dz);

		while(node >= 0)
    {
				if(node < All.MaxPart)  // single particle
				{
						p = node;
						node = Nextnode[node];

						if(P[p].Type > 0) // not gas particle
							continue;

						//if(P[p].Ti_current != All.Ti_Current)
						//	drift_particle(p, All.Ti_Current);

						dx = P[p].Pos[0] - pt.x;
						if(dx > cur_mindist)
								continue;
								
						dy = P[p].Pos[1] - pt.y;
						if(dy > cur_mindist)
								continue;
							
						dz = P[p].Pos[2] - pt.z;
						if(dz > cur_mindist)
								continue;
							
						double curdist2 = dx * dx + dy * dy + dz * dz;
						if(curdist2 > cur_mindist * cur_mindist)
								continue;

						cur_mindist = sqrt(curdist2);
						nearest = p;
				}
				else
				{
						if(node >= All.MaxPart + MaxNodes) { // pseudo particle
								// what's going on here is fuzzy to me
								node = Nextnode[node - MaxNodes];
								continue;
						}

						current = &Nodes[node];

						//if(current->Ti_current != All.Ti_Current)
						//force_drift_node(node, All.Ti_Current);

						if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
						{
								if(current->u.d.mass) // open cell
								{
										node = current->u.d.nextnode;
										continue;
								}
						}

						// in case the node can be discarded
						node = current->u.d.sibling;

						// first quick tests along the axes
						double test_dist = cur_mindist + 0.5 * current->len;
						
						dx = current->center[0] - pt.x;
						if(dx > test_dist)
								continue;
								
						dy = current->center[1] - pt.y;
						if(dy > test_dist)
								continue;
								
						dz = current->center[2] - pt.z;
						if(dz > test_dist)
								continue;

						// now test against the minimal sphere enclosing everything
						test_dist += FACT1 * current->len;
						if(dx * dx + dy * dy + dz * dz > test_dist * test_dist)
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

bool ArepoMesh::AdvanceRayOneCell(const Ray &ray, float *t0, float *t1, Spectrum &Lv, Spectrum &Tr)
{
		int qmin = -1;
		
		if (ray.task != ThisTask) {
				cout << "[" << ThisTask << "] ERROR! ray.task = " << ray.task << " differs from ThisTask!" << endl;
				endrun(1117);
		}
		
		int q        = SphP[ray.index].first_connection;
		double min_t = MAX_REAL_NUMBER;

		int dp;
		Vector midp, norm;		
		
		IF_DEBUG(cout << " starting face q = " << q << endl);
		
		Point hitbox  = ray(*t0);
		Point exitbox = ray(*t1);
		
		// examine faces to find exit point
		while (q >= 0)
		{
				dp = DC[q].dp_index;
				
				// midpoint
				midp = Vector( 0.5 * (DP[dp].x + P[ray.index].Pos[0]),
				               0.5 * (DP[dp].y + P[ray.index].Pos[1]),
									     0.5 * (DP[dp].z + P[ray.index].Pos[2]) );
				
				// vector pointing to the outside, normal to a voronoi face of the cell
				norm = Vector( (DP[dp].x - P[ray.index].Pos[0]),
											 (DP[dp].y - P[ray.index].Pos[1]),
											 (DP[dp].z - P[ray.index].Pos[2]) );
				
				//IF_DEBUG(midp.print(" midp "));
				//IF_DEBUG(norm.print(" norm "));
				
				// find intersection of ray with this face
				//double dotprod = Dot( (exitbox - hitbox), norm );
				double dotprod = Dot( ray.d, norm );
				
				if (dotprod > 0) {
						//dist of intersection from ray.o (any point on the ray)
						double t = ((midp.x - hitbox.x) * norm.x +
											  (midp.y - hitbox.y) * norm.y +
											  (midp.z - hitbox.z) * norm.z)	/ dotprod;
						//IF_DEBUG(cout << " q[" << q << "] dotprod>0 and t = " << t << " (min=" << (ray.min_t-*t0) << ")" << endl);
						
						if (t > (ray.min_t-*t0) && t < min_t) {
								IF_DEBUG(cout << " intersection t = " << t << " setting new min_t, qmin = q = " << q << endl);
								min_t = t;
								qmin = q;
						}
				}
		
				// move to next face
				if (q == SphP[ray.index].last_connection) {
						//IF_DEBUG(cout << " q = last_connection = " << q << " (done with faces)" << endl);
						break;
				}
						
				//IF_DEBUG(cout << " done with q = " << q << " (dp=" << dp << ") next = " << DC[q].next << endl);
				q = DC[q].next;
		}
		
		// check for proper exit point
		if (qmin != -1) {
		
				// clamp min_t to avoid integrating outside the box
				IF_DEBUG(cout << " min_t = " << min_t << " (t1=" << *t1 << " t0=" << *t0 << ")" << endl);
				min_t = Clamp(min_t,0.0,(*t1-*t0));
				
				Point hitcell  = ray(ray.min_t);
				Point exitcell = ray(*t0 + min_t);
				
				IF_DEBUG(hitcell.print(" hcell "));
				IF_DEBUG(exitcell.print(" ecell "));

				// find interpolated value at midpoint of line segment through voronoi cell
				// or sample at stepsize through cell (needed for smooth transfer functions)
		
				Vector sphCen = Vector(SphP[ray.index].Center[0],
															 SphP[ray.index].Center[1],
															 SphP[ray.index].Center[2]);
				Vector sphDensGrad = Vector(SphP[ray.index].Grad.drho[0],
																		SphP[ray.index].Grad.drho[1],
																		SphP[ray.index].Grad.drho[2]);
														 
				norm = exitcell - hitcell;

				IF_DEBUG(norm.print(" norm "));
				IF_DEBUG(sphCen.print(" sphCen "));
				
				// compute total path length through cell
				double len = norm.Length();
				
				// TODO: sub-step length should be adaptive based on gradients
				//       should only substep if the TF will evaluate nonzero somewhere inside
				
				// decide on sub-stepping
				if (viStepSize && len > viStepSize) {
						// sub-stepping (cell too large), get a number of values along the segment
						int nSamples = (int)ceilf(len / viStepSize);
						double fracstep = 1.0 / nSamples;
						double halfstep = 0.5 / nSamples;
						
						IF_DEBUG(cout << " sub-stepping len = " << len << " nSamples = " << nSamples 
													<< " (step = " << len/nSamples << ")" << endl);
						
						for (int i = 0; i < nSamples; ++i) {
								Vector midpt = Vector(hitcell[0] + (i*fracstep+halfstep) * norm[0],
																			hitcell[1] + (i*fracstep+halfstep) * norm[1],
																			hitcell[2] + (i*fracstep+halfstep) * norm[2]) - sphCen;
								
								IF_DEBUG(midpt.print(" substep midpt rel sphcen "));
								
								double rho = SphP[ray.index].Density + Dot(sphDensGrad,midpt);
								double utherm = 0.0;
								
								// multiplying by the stepsize makes these integrated (total) quantities, e.g. for 
								// a "projected column density" image. for sampling quantities with a transfer 
								// function we just want their true value at every point in space
								// TODO: make this multiplication optional based on Config
								//rho    *= step;
								//utherm *= step;
								
								IF_DEBUG(cout << "  segment[" << i << "] fractrange [" << (i*fracstep) << "," 
															<< (i*fracstep)+fracstep << "] rho = " << SphP[ray.index].Density
															<< " rho w/ grad = " << rho << " (grad.x = " << sphDensGrad.x
															<< " grad.y = " << sphDensGrad.y << " grad.z = " << sphDensGrad.z << ")" << endl);
								
								// TODO tau calculation
								
								// Compute emission-only source term with transfer function
								float vals[2]; vals[0] = (float)rho; vals[1] = (float)utherm;
								
								Lv += Tr * transferFunction->Lve(vals);
						}
						
				} else {
						// no sub-stepping (cell sufficiently small), get interpolated values at segment center
						Vector midpt = Vector(hitcell[0] + 0.5 * norm[0],
																	hitcell[1] + 0.5 * norm[1],
																	hitcell[2] + 0.5 * norm[2]) - sphCen;
						
						IF_DEBUG(midpt.print(" onestep midp rel sphcen "));
						
						double rho = SphP[ray.index].Density + Dot(sphDensGrad,midpt);
						double utherm = 0.0;
						
						//rho    *= len;
						//utherm *= len;
						
						IF_DEBUG(cout << " segment len = " << len << " rho = " << SphP[ray.index].Density 
													<< " grad.x = " << sphDensGrad.x << " rho w/ grad = " << rho << endl);

						//TODO tau calculation
						//Spectrum stepTau = scene->volumeRegion->tau(tauRay,0.5f * stepSize, rng.RandomFloat());
						//Tr *= Exp(-stepTau);

						// Compute emission-only source term with transfer function
						float vals[2]; vals[0] = (float)rho; vals[1] = (float)utherm;
						
						Lv += Tr * transferFunction->Lve(vals);
				}
				
				// update ray: transfer to next voronoi cell (possibly on different task)
				ray.task  = DC[qmin].task;
				ray.index = DC[qmin].index;
				ray.min_t = Clamp(min_t + *t0,ray.min_t,ray.max_t);
				
				IF_DEBUG(cout << " updated ray new task = " << ray.task << " index = " << ray.index 
											<< " (" << DC[qmin].index << ") min_t = " << ray.min_t << endl);
				
				if (fabs(ray.min_t - ray.max_t) <= INSIDE_EPS) {
						// apparently this ray is done?
						IF_DEBUG(cout << " min_t == t1 = " << *t1 << ", ray done." << endl);
						return false;
				}
				if (ray.index == -1) {
						// no connection on other side of face?
						cout << " ERROR. next cell index = -1" << " (ray.min_t=" << ray.min_t 
								 << " max_t=" << ray.max_t << " diff = " << fabs(ray.min_t - ray.max_t) << ")" << endl;
						endrun(1109);
				}
		} else {
				// apparently this ray is done?
				cout << " failed to intersect a face, hope this ray is done?" << endl;
				return false;
		}

		return true;
}

void ArepoMesh::ComputeQuantityBounds()
{
		float pmax  = 0.0;
		float pmin  = INFINITY;
		float pmean = 0.0;
		
		for (int i=0; i < N_gas; i++) {
				if (SphP[i].Density > pmax)
						pmax = SphP[i].Density;
				if (SphP[i].Density < pmin)
						pmin = SphP[i].Density;
				pmean += SphP[i].Density;
		}
		
		densBounds[0] = pmin;
		densBounds[1] = pmax;
		densBounds[2] = pmean / N_gas;
		
		IF_DEBUG(cout << " Density min = " << densBounds[0] << " max = " << densBounds[1] 
									<< " mean = " << densBounds[2] << endl);

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
