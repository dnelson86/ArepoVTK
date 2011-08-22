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
}

void Arepo::Cleanup()
{
		MPI_Finalize();
}

bool Arepo::LoadSnapshot()
{
    IF_DEBUG(cout << "Arepo::LoadSnapshot(" << snapFilename << ")." << endl);
		
		freopen("/dev/null","w",stdout); //hide arepo stdout
		freopen("/dev/null","w",stderr); //hide MPI errors (TODO:temp)
		
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
		freopen("/dev/tty","w",stderr);
		
		return true;
}

ArepoMesh::ArepoMesh(const TransferFunction *tf, const Transform &VolumeToWorld)
{
		IF_DEBUG(cout << "ArepoMesh() constructor." << endl);
		
		// transfer function and transform
		transferFunction = tf;
		viStepSize       = Config.viStepSize;
		WorldToVolume    = Inverse(VolumeToWorld);
		
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
		
		// boxsize
		extent = BBox(Point(0.0,0.0,0.0),Point(boxSize,boxSize,boxSize));
		
		// debugging
		//arepoMesh->DumpMesh();
}

void ArepoMesh::LocateEntryCell(const Ray &ray, float *t0, float *t1)
{
		double pos[3], mindist;
		Point hitbox  = WorldToVolume(ray(*t0));
		Point exitbox = WorldToVolume(ray(*t1));
		
		pos[0] = hitbox.x;
		pos[1] = hitbox.y;
		pos[2] = hitbox.z;
		
		IF_DEBUG(cout << " ray hits box at W x = " << ray(*t0).x << " y = " << ray(*t0).y << " z = " << ray(*t0).z << endl);
		IF_DEBUG(cout << " ray hits box at V x = " << hitbox.x << " y = " << hitbox.y << " z = " << hitbox.z << endl);
		IF_DEBUG(cout << " ray exit box at W x = " << ray(*t1).x << " y = " << ray(*t1).y << " z = " << ray(*t1).z << endl);
		IF_DEBUG(cout << " ray exit box at V x = " << exitbox.x << " y = " << exitbox.y << " z = " << exitbox.z << endl);
		
		// use tree-finding function to find closest gas point (local only)
		const int corig2 = ngb_treefind_nearest_local(&pos[0], &mindist); //,0
		
		int corig;
		mindist = MAX_REAL_NUMBER;
		for (int i=0; i < Ndp; i++) {
				double dist = sqrt( (pos[0]-DP[i].x)*(pos[0]-DP[i].x) + 
													  (pos[1]-DP[i].y)*(pos[1]-DP[i].y) + 
														(pos[2]-DP[i].z)*(pos[2]-DP[i].z) );
				//cout << " DP[" << i << "] dist = " << dist << endl;
				if (dist < mindist) {
						mindist = dist;
						corig = i;
				}
		}
		
		IF_DEBUG(cout << " corig = " << corig << " DP[corig].index = " << DP[corig].index 
									<< " dist = " << mindist << " (x = " << DP[corig].x << " y = "
									<< DP[corig].y << " z = " << DP[corig].z << ")" << endl);
		IF_DEBUG(cout << " corig2 = " << corig2 << " DP[corig2].index = " << DP[corig2].index 
									<< " dist = " << mindist << " (x = " << DP[corig2].x << " y = "
									<< DP[corig2].y << " z = " << DP[corig2].z << ")" << endl);						
				 
		// refine nearest point search
/*		int cc, c = corig;
		
		for (i = 0; i <= 4; i++) {
				cc = find_closest_neighbor_in_cell(p,c); //TODO
				if (cc = c)
						break;
				c = cc;
		} */
		
		ray.index = corig; //cc
		//cout << " iterates i = " << i << " picked ray.index = " << ray.index << endl;
		
		if (ray.index < 0 || ray.index >= N_gas) {
				cout << "ERROR! ray.index = " << ray.index << " out of bounds." << endl;
				endrun(1107);
		}
		
}

bool ArepoMesh::AdvanceRayOneCell(const Ray &ray, float *t0, float *t1, Spectrum &Lv, Spectrum &Tr)
{
		int qmin = -1;
		
		int q        = SphP[ray.index].first_connection;
		double min_t = MAX_REAL_NUMBER;

		int dp;
		Vector midp, norm;		
		
		IF_DEBUG(cout << " starting face q = " << q << endl);
		
		Point hitbox  = WorldToVolume(ray(*t0));
		Point exitbox = WorldToVolume(ray(*t1));		
		
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
				
				IF_DEBUG(midp.print(" midp "));
				IF_DEBUG(norm.print(" norm "));
				
				// find intersection of ray with this face
				//double dotprod = Dot( (exitbox - hitbox), norm );
				double dotprod = Dot( ray.d, norm );
				
				if (dotprod > 0) {
						//dist of intersection from ray.o (any point on the ray)
						double t = ((midp.x - hitbox.x) * norm.x +
											  (midp.y - hitbox.y) * norm.y +
											  (midp.z - hitbox.z) * norm.z)	/ dotprod;
						//double t = Dot( Vector(midp.x - ray.o.x,
						//                       midp.y - ray.o.y,
						//											 midp.z - ray.o.z), norm ) / dotprod;
						IF_DEBUG(cout << " q[" << q << "] dotprod>0 and t = " << t << " (min=" << (ray.min_t-*t0) << ")" << endl);
						
						if (t > (ray.min_t-*t0) && t < min_t) {
								IF_DEBUG(cout << " intersection V t = " << t << " setting new min_t, qmin = q = " << q << endl);
								min_t = t; //V
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
				min_t = Clamp(min_t,0.0,(*t1-*t0)); //V
				
				IF_DEBUG(ray(ray.min_t-*t0).print(" hcell W "));
				IF_DEBUG(ray(min_t).print(" ecell W "));
				
				IF_DEBUG(WorldToVolume(ray(ray.min_t)).print(" hcell V "));
				IF_DEBUG(WorldToVolume(ray(min_t+*t0)).print(" ecell V "));

				// find interpolated value at midpoint of line segment through voronoi cell
				// or sample at stepsize through cell (needed for smooth transfer functions)
		
				Vector sphCen = Vector(SphP[ray.index].Center[0],
															 SphP[ray.index].Center[1],
															 SphP[ray.index].Center[2]);
				Vector sphDensGrad = Vector(SphP[ray.index].Grad.drho[0],
																		SphP[ray.index].Grad.drho[1],
																		SphP[ray.index].Grad.drho[2]);
														 
				norm = exitbox - hitbox;

				IF_DEBUG(norm.print(" norm V "));
				IF_DEBUG(sphCen.print(" sphCen V "));
				IF_DEBUG(midp.print(" midp V "));
				
				// compute total path length through cell
				double len = (min_t - (ray.min_t-*t0)) * norm.Length();
				
				// decide on sub-stepping
				if (viStepSize && len > viStepSize) {
						// sub-stepping (cell too large), get a number of values along the segment
						int nSamples = (int)ceilf(len / viStepSize);
						double step = len / nSamples;
						double t0step = *t0;
						
						IF_DEBUG(cout << " sub-stepping len = " << len << " nSamples = " << nSamples 
													<< " (step = " << step << ")" << endl);
						
						for (int i = 0; i < nSamples; ++i, t0step += step) {
								midp = Vector(hitbox) + 0.5 * (ray.min_t-*t0 + t0step) * norm - sphCen;
								
								double rho = SphP[ray.index].Density + Dot(sphDensGrad,midp);
								double utherm = 0.0;
								
								// multiplying by the stepsize makes these integrated (total) quantities, e.g. for 
								// a "projected column density" image. for sampling quantities with a transfer 
								// function we just want their true value at every point in space
								// TODO: make this multiplication optional based on Config
								//rho    *= step;
								//utherm *= step;
								
								IF_DEBUG(cout << "  segment[" << i << "] trange [" << (ray.min_t-*t0 + t0step) << "," 
															<< (ray.min_t-*t0 + t0step)+step << "] rho = " << SphP[ray.index].Density*step
															<< " rho w/ grad = " << rho << " (grad.x = " << sphDensGrad.x
															<< " grad.y = " << sphDensGrad.y << " grad.z = " << sphDensGrad.z << ")" << endl);
								
								// TODO tau calculation
								
								// Compute emission-only source term with transfer function
								float vals[2]; vals[0] = (float)rho; vals[1] = (float)utherm;
								
								Lv += Tr * transferFunction->Lve(vals);
						}
						
				} else {
						// no sub-stepping (cell sufficiently small), get interpolated values at segment center
						midp = Vector(hitbox) + 0.5 * (ray.min_t-*t0 + min_t) * norm - sphCen;
						//midp = Vector(hitbox) + 0.5 * (min_t - *t0) * norm - sphCen;
						
						double rho = SphP[ray.index].Density + Dot(sphDensGrad,midp);
						double utherm = 0.0;
						
						// multiplying by the stepsize makes these integrated (total) quantities, e.g. for 
						// a "projected column density" image. for sampling quantities with a transfer 
						// function we just want their true value at every point in space
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
				ray.min_t = Clamp(min_t + *t0,ray.min_t,ray.max_t); //W [t0,t1]
				
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
					
				pts[j] = Inverse(WorldToVolume)(Point(DP[DT[i].p[j]].x,
																						  DP[DT[i].p[j]].y, 
																						  DP[DT[i].p[j]].z));
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
				
				edges->push_back(Line(Inverse(WorldToVolume)(prev),Inverse(WorldToVolume)(next)));
				
				IF_DEBUG(cout << " edge[" << i_face << "," << i 
											<< "] prev.x = " << prev.x << " prev.y = " << prev.y << " prev.z = " << prev.z
											<<  " next.x = " << next.x << " next.y = " << next.y << " next.z = " << next.z << endl);
						 
				prev = next;
		}
		
		return true;
}

ArepoMesh *CreateArepoMesh(const TransferFunction *tf, const Transform &volume2world)
{
    return new ArepoMesh(tf, volume2world);
}

#endif //ENABLE_AREPO
