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
		
#ifndef DEBUG
		freopen("/dev/null","w",stdout); //hide arepo stdout
#endif		
		
		// set startup options
		WriteMiscFiles = 0;
		RestartSnapNum = -1;
		RestartFlag    = SUNRISE;

		strcpy(ParameterFile,paramFilename.c_str());

		// call arepo: run setup
		begrun1();
		
		// check snapshot exists
		if (!ifstream(snapFilename.c_str())) {
				freopen("/dev/tty","w",stdout);
				cout << "Arepo::LoadSnapshot() ERROR! Snapshot [" << snapFilename << "] not found." << endl;
				endrun(1121);
		}
		
		if (snapFilename.substr(snapFilename.size()-5,5) == ".hdf5") // cut off extension
		  snapFilename = snapFilename.substr(0,snapFilename.size()-5);

		// load snapshot
		read_ic(snapFilename.c_str());
		
		// call arepo: read snapshot, allocate storage for tree, 
		//             initialize particle data, domain decomposition, initial HSML
  	if (init() != SUNRISE) {
				cout << "Arepo::LoadSnapshot() ERROR: Arepo did not return successfully." << endl;
				return false;
		}

#ifndef DEBUG
		//freopen("/dev/tty","w",stdout); //return stdout to terminal
		freopen("LM.out","w",stdout); //return stdout to a file (LSF) //TODO
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
		
		// preprocessing
		//ArepoMesh::ComputeVoronoiEdges();
		ArepoMesh::ComputeQuantityBounds();
		ArepoMesh::CalculateMidpoints();
		
		// boxsize
		extent = BBox(Point(0.0,0.0,0.0),Point(boxSize,boxSize,boxSize));
		
		IF_DEBUG(extent.print(" ArepoMesh extent "));
		
		// TODO: temp units
		unitConversions[TF_VAL_DENS]   = All.UnitDensity_in_cgs / MSUN_PER_PC3_IN_CGS;
		unitConversions[TF_VAL_UTHERM] = All.UnitEnergy_in_cgs;
		
		IF_DEBUG(cout << "unitConv[dens]   = " << unitConversions[TF_VAL_DENS] << endl);
		IF_DEBUG(cout << "unitConv[utherm] = " << unitConversions[TF_VAL_UTHERM] << endl);
		
		// debugging
		ArepoMesh::DumpMesh();
}

void ArepoMesh::LocateEntryCellBrute(const Ray &ray, float *t0, float *t1)
{
		// note: using the brute force search is O(N_rays * NumPart) - not good
		Point hitbox  = ray(*t0);
		
		int corig;
		int corig2;
		double mindist = MAX_REAL_NUMBER;
		double mindist2 = MAX_REAL_NUMBER;
		
		for (int i=0; i < N_gas; i++) {
				double dist = sqrt( (hitbox.x-P[i].Pos[0])*(hitbox.x-P[i].Pos[0]) + 
													  (hitbox.y-P[i].Pos[1])*(hitbox.y-P[i].Pos[1]) + 
														(hitbox.z-P[i].Pos[2])*(hitbox.z-P[i].Pos[2]) );
				if (dist < mindist) {
						mindist = dist;
						corig = i;
				}
		}
		
		for (int i=0; i < Ndp; i++) {
				double dist = sqrt( (hitbox.x-DP[i].x)*(hitbox.x-DP[i].x) + 
													  (hitbox.y-DP[i].y)*(hitbox.y-DP[i].y) + 
														(hitbox.z-DP[i].z)*(hitbox.z-DP[i].z) );
				if (dist < mindist2) {
						mindist2 = dist;
						corig2 = i;
				}
		}

		IF_DEBUG(cout << " brute SphP corig = " << corig << " P[corig].index = none" 
									<< " dist = " << mindist << " (x = " << P[corig].Pos[0] << " y = "
									<< P[corig].Pos[1] << " z = " << P[corig].Pos[2] << ")" << endl);
		
		if (corig2 >= N_gas) {
				cout << " BRUTE DP IS LOCAL GHOST corig2 = " << corig2 << " (N_gas = " << N_gas << ")" << endl;
		}
		
		IF_DEBUG(cout << " brute DP corig = " << corig2 << " P[corig].index = none" 
									<< " dist = " << mindist2 << " (x = " << P[corig2].Pos[0] << " y = "
									<< P[corig2].Pos[1] << " z = " << P[corig2].Pos[2] << ")" << endl);
		
		if (corig != corig2) {
				cout << " WARNING: BRUTE SPHP != DP" << endl;
				endrun(1135);
		}
		
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
		
		// check indexing
		if (DP[ray.index].index != ray.index) {
				cout << "ERROR! Chosen tetra does not link back to nearest gas particle " 
						 << ray.index << " " << DP[ray.index].index << endl;
				endrun(1129);
		}
							
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
								node = Nextnode[node - MaxNodes];
								continue;
						}

						current = &Nodes[node];

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
		
		// check for local ghosts
		if (ray.index >= N_gas) {
				cout << "ERROR! ray.index = " << ray.index << " (N_gas=" << N_gas << " Ndp=" << Ndp 
						 << ") [task=" << ray.task << "] out of bounds." << endl;
				cout << " hit pos: " << hitbox.x << " " << hitbox.y << " " << hitbox.z << endl;
				cout << " DP[ray.index] pos: " << DP[ray.index].x << " " << DP[ray.index].y << " " << DP[ray.index].z << endl;
				cout << " P[ray.index-N_gas] (local) pos: " << P[ray.index-N_gas].Pos[0] << " " << P[ray.index-N_gas].Pos[1]
						 << " " << P[ray.index-N_gas].Pos[2] << endl;
				endrun(1107);
		}
		
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
				
				// find intersection of ray with this face
				double dotprod = Dot( ray.d, norm );
				
				if (dotprod > 0) {
						//dist of intersection from ray.o (any point on the ray)
						double t = ((midp.x - hitbox.x) * norm.x +
											  (midp.y - hitbox.y) * norm.y +
											  (midp.z - hitbox.z) * norm.z)	/ dotprod;
						IF_DEBUG(cout << " q[" << q << "] dotprod>0 and t = " << t << " (min=" << (ray.min_t-*t0) << ")" << endl);
						
						if (t > (ray.min_t-*t0) && t < min_t) {
								IF_DEBUG(cout << " intersection t = " << t << " setting new min_t, qmin = q = " << q << endl);
								min_t = t;
								qmin = q;
						}
				}
		
				// move to next face
				if (q == SphP[ray.index].last_connection) {
						IF_DEBUG(cout << " q = last_connection = " << q << " (done with faces)" << endl);
						break;
				}
						
				IF_DEBUG(cout << " done with q = " << q << " (dp=" << dp << ") next = " << DC[q].next << endl);
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

				// cell gradient
				Vector sphCen(SphP[ray.index].Center[0],
											SphP[ray.index].Center[1],
											SphP[ray.index].Center[2]);
				Vector sphDensGrad(SphP[ray.index].Grad.drho[0],
													 SphP[ray.index].Grad.drho[1],
													 SphP[ray.index].Grad.drho[2]);
														 
				norm = exitcell - hitcell;

				IF_DEBUG(norm.print(" norm "));
				IF_DEBUG(sphCen.print(" sphCen "));
				
				// compute total path length through cell
				double len = norm.Length();
				
				// find interpolated density at midpoint of line segment through voronoi cell
				Vector midpt(hitcell[0] + 0.5 * norm[0],
										 hitcell[1] + 0.5 * norm[1],
										 hitcell[2] + 0.5 * norm[2]);
        midpt -= sphCen;
				
				IF_DEBUG(midpt.print(" onestep midp rel sphcen "));
				
				double rho    = SphP[ray.index].Density;
				double utherm = SphP[ray.index].Utherm;
				
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
				Tr *= Exp(-stepTau);
				
				// TODO: sub-step length should be adaptive based on gradients
				// TODO: should only substep if the TF will evaluate nonzero somewhere inside
				
				// decide on sub-stepping
				if (viStepSize && len > viStepSize) {
						// sub-stepping (cell too large), get a number of values along the segment
						int nSamples = (int)ceilf(len / viStepSize);
						double fracstep = 1.0 / nSamples;
						double halfstep = 0.5 / nSamples;
						
						IF_DEBUG(cout << " sub-stepping len = " << len << " nSamples = " << nSamples 
													<< " (step = " << len/nSamples << ")" << endl);
						
						for (int i = 0; i < nSamples; ++i) {
								Vector midpt(hitcell[0] + (i*fracstep+halfstep) * norm[0],
														 hitcell[1] + (i*fracstep+halfstep) * norm[1],
														 hitcell[2] + (i*fracstep+halfstep) * norm[2]);
                midpt -= sphCen;
								
								IF_DEBUG(midpt.print(" substep midpt rel sphcen "));
								
								double rho    = SphP[ray.index].Density + Dot(sphDensGrad,midpt);
								
								IF_DEBUG(cout << "  segment[" << i << "] fractrange [" << (i*fracstep) << "," 
															<< (i*fracstep)+fracstep << "] rho = " << SphP[ray.index].Density
															<< " rho w/ grad = " << rho << " (grad.x = " << sphDensGrad.x
															<< " grad.y = " << sphDensGrad.y << " grad.z = " << sphDensGrad.z << ")" << endl);
								
								// pack cell quantities for TF
								float vals[TF_NUM_VALS];
								
								vals[TF_VAL_DENS]        = (float) rho;
								vals[TF_VAL_UTHERM]      = (float) SphP[ray.index].Utherm; //TODO: gradient (MATERIALS)
								vals[TF_VAL_PRES]        = (float) SphP[ray.index].Pressure; //TODO: gradient
								vals[TF_VAL_ENERGY]      = (float) SphP[ray.index].Energy;
								
								vals[TF_VAL_VEL_X]       = (float) P[ray.index].Vel[0]; //TODO: gradients
								vals[TF_VAL_VEL_Y]       = (float) P[ray.index].Vel[1];
								vals[TF_VAL_VEL_Z]       = (float) P[ray.index].Vel[2];
								vals[TF_VAL_VEL_DIV]     = (float) SphP[ray.index].DivVel;
								vals[TF_VAL_VEL_CURL]    = (float) SphP[ray.index].CurlVel;
								
								vals[TF_VAL_POTENTIAL]   = (float) P[ray.index].Potential;
								
#ifdef METALS
								vals[TF_VAL_METALLICITY] = (float) SphP[ray.index].Metallicity;
#endif
#ifdef COOLING
								vals[TF_VAL_NE]          = (float) SphP[ray.index].Ne;
#endif
#ifdef USE_SFR
								vals[TF_VAL_SFR]         = (float) SphP[ray.index].Sfr;
#endif
								
								// compute emission-only source term using transfer function
								Lv += Tr * transferFunction->Lve(vals);
						}
						
				} else {
						// no sub-stepping (cell sufficiently small), get interpolated values at segment center

						// apply TF to integrated (total) quantities (only appropriate for constant?)
						if (Config.projColDens) {
								rho    *= len;
								utherm *= len;
						}
						
						IF_DEBUG(cout << " segment len = " << len << " rho = " << SphP[ray.index].Density 
													<< " grad.x = " << sphDensGrad.x << " rho w/ grad = " << rho << endl);

						// pack cell quantities for TF
						float vals[TF_NUM_VALS];
						
								vals[TF_VAL_DENS]        = (float) rho;
								vals[TF_VAL_UTHERM]      = (float) utherm; //TODO: gradient (MATERIALS)
								vals[TF_VAL_PRES]        = (float) SphP[ray.index].Pressure; //TODO: gradient
								vals[TF_VAL_ENERGY]      = (float) SphP[ray.index].Energy;
								
								vals[TF_VAL_VEL_X]       = (float) P[ray.index].Vel[0]; //TODO: gradients
								vals[TF_VAL_VEL_Y]       = (float) P[ray.index].Vel[1];
								vals[TF_VAL_VEL_Z]       = (float) P[ray.index].Vel[2];
								vals[TF_VAL_VEL_DIV]     = (float) SphP[ray.index].DivVel;
								vals[TF_VAL_VEL_CURL]    = (float) SphP[ray.index].CurlVel;
								
								vals[TF_VAL_POTENTIAL]   = (float) P[ray.index].Potential;
								
#ifdef METALS
								vals[TF_VAL_METALLICITY] = (float) SphP[ray.index].Metallicity;
#endif
#ifdef COOLING
								vals[TF_VAL_NE]          = (float) SphP[ray.index].Ne;
#endif
#ifdef USE_SFR
								vals[TF_VAL_SFR]         = (float) SphP[ray.index].Sfr;
#endif
						
						// compute emission-only source term using transfer function
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
				// failed to intersect a face (shouldn't really happen?)
#ifdef DEBUG
				cout << " failed to intersect a face, hope this ray is done?" << endl;
				if (ray.min_t < ray.max_t - INSIDE_EPS) {
						cout << "ERROR! Ray did not finish. min_t = " << ray.min_t << " max_t = " << ray.max_t << endl;
						cout << " ray pos: " << hitbox.x << " " << hitbox.y << " " << hitbox.z << endl;
						cout << " DP[ray.index] pos: " << DP[ray.index].x << " " << DP[ray.index].y 
								 << " " << DP[ray.index].z << endl << endl;
						//endrun(1130);
				}
#endif
				return false;
		}

		return true;
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
						
						if (DP[i].index >= 0 && DP[i].index < N_gas)
								SphP_ID_n = DP[i].index;
						else if (i >= N_gas)
								SphP_ID_n = DP[i].index - N_gas;
								
						// skip invalid neighbors
						if (SphP_ID_n < 0)
								continue;
								
						const Vector midp( 0.5 * (cellp + Vector(DP[SphP_ID_n].x,
																										 DP[SphP_ID_n].y,
																										 DP[SphP_ID_n].z)) );
																										 
						// add connection to the connectivity map
						midpoints.push_back(midp);
						opposite_points.push_back(SphP_ID_n);
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
