/*
 * arepoInterp.cpp
 * dnelson
 */
 
#include <alloca.h>

#include "transform.h"
#include "spectrum.h"
#include "arepo.h"
#include "volume.h"
#include "transfer.h"

#ifdef ENABLE_AREPO

// NNI_WATSON_SAMBRIDGE
#define MAX_NUM_TETS 100
#define MAX_NUM_NODES 400
#define NODE_A 3
#define NODE_B 0
#define NODE_C 1
#define NODE_D 2

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
						
				//cout << "[" << i << "] p0 = " << DT[i].p[0] << " p1 = " << DT[i].p[2] << 
				//	" p2 = " << DT[i].p[2] << " p3 = " << DT[i].p[3] << " det_A = " << det_A << endl;
		
				if (det_A < INSIDE_EPS && det_A > -INSIDE_EPS) { // debugging
				  //cout << "ERROR: det_A is zero! " << det_A << " [" << i << "]" << endl;
					//exit(12509);
					continue;
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

float ArepoMesh::calcNeighborHSML(int sphInd, Vector &pt, int dpInd)
{
		float dx,dy,dz,xtmp,ytmp,ztmp;
		float distsq, hsml2 = 0.0;
		
#ifdef USE_DC_CONNECTIVITY
		
		int edge = SphP[sphInd].first_connection;
		int last_edge = SphP[sphInd].last_connection;
		
		while(edge >= 0) {
		  int sphp_neighbor = DC[edge].index;
			
			// could connect to bounding tetra if we don't have a ghost across this boundary
			if ( sphp_neighbor < 0
#ifdef NO_GHOST_CONTRIBS
			|| DC[edge].dp_index >= NumGas
#endif
			) {
				if ( edge == last_edge )
					break;
				edge = DC[edge].next;
				continue;
			}
			
			dx = NGB_PERIODIC_LONG_X(P[sphp_neighbor].Pos[0] - pt.x);
			dy = NGB_PERIODIC_LONG_Y(P[sphp_neighbor].Pos[1] - pt.y);
			dz = NGB_PERIODIC_LONG_Z(P[sphp_neighbor].Pos[2] - pt.z);
			distsq = dx*dx + dy*dy + dz*dz;
			  
      if(distsq > hsml2)
        hsml2 = distsq;
				
			// move to next neighbor
			if(edge == last_edge)
				break;
				
#ifdef DEBUG
			if (DC[edge].next == edge || DC[edge].next < 0)
			  terminate(" what is going on (%d %d %d %d)",DC[edge].next,edge,DC[edge].dp_index,sphp_neighbor);
#endif				

			edge = DC[edge].next;
		}
		
#else // USE_ALTERNATIVE_CONNECTIVITY
		
		const int start_edge = midpoint_idx[dpInd].first;
		const int n_edges    = midpoint_idx[dpInd].second;
		
		for ( int k=0; k < n_edges; k++ ) {
			int sphp_neighbor = getSphPID( opposite_points[start_edge+k] );
			
			dx = NGB_PERIODIC_LONG_X(P[sphp_neighbor].Pos[0] - pt.x);
			dy = NGB_PERIODIC_LONG_Y(P[sphp_neighbor].Pos[1] - pt.y);
			dz = NGB_PERIODIC_LONG_Z(P[sphp_neighbor].Pos[2] - pt.z);
			distsq = dx*dx + dy*dy + dz*dz;
			  
      if(distsq > hsml2)
        hsml2 = distsq;
		}

#endif // CONNECTIVITY

    float hinv = 1.0 / sqrtf(hsml2);
		
		return hinv;

}
#endif

#ifdef NATURAL_NEIGHBOR_INTERP
void inline periodic_wrap_DP_point(point &dp_pt, Vector &ref)
{
#ifdef PERIODIC
  double dx, dy, dz;

  dx = dp_pt.x - ref.x;
  dy = dp_pt.y - ref.y;
  dz = dp_pt.z - ref.z;

  if(dx > boxHalf_X)
    dp_pt.x -= boxSize_X;
  if(dx < -boxHalf_X)
    dp_pt.x += boxSize_X;
  if(dy > boxHalf_Y)
    dp_pt.y -= boxSize_Y;
  if(dy < -boxHalf_Y)
    dp_pt.y += boxSize_Y;
  if(dz > boxHalf_Z)
    dp_pt.z -= boxSize_Z;
  if(dz < -boxHalf_Z)
    dp_pt.z += boxSize_Z;
#endif
}
#endif

#ifdef NNI_WATSON_SAMBRIDGE
inline bool ArepoMesh::needTet(int tt, point *pp, int *tet_inds, int *nTet)
{
	// check if this tetra is already in tet_inds, if so skip
	// we do this instead of a "mark" (could use DT[tt].s[0]) to make this thread safe
	for (int i=0; i < *nTet; i++)
	  if ( tet_inds[i] == tt )
		  return false;

	// check if pp is inside tetra tt, use non-exact if possible
	int ret;
	
	ret = InSphere_Errorbound( &DP[DT[tt].p[0]], &DP[DT[tt].p[1]], &DP[DT[tt].p[2]],
														 &DP[DT[tt].p[3]], pp ); // -1 outside, +1 inside, 0 need exact

	//IF_DEBUG(cout << "   needTet [" << tt << "] ret = " << ret << " pp.x = " << pp->x << endl);
	
	if ( ret == 0 ) {
		terminate(" handle the request for exact insphere ");
	}
	else if ( ret > 0 ) {
		return true; // add
	}
	
	return false; // skip
}

void ArepoMesh::addTet(int tt, point *pp, int *node_inds, int *tet_inds, int *nNode, int *nTet)
{
	if ( *nNode+3 >= MAX_NUM_TETS || *nTet >= MAX_NUM_NODES )
	  terminate(" addTet reached maximum buffer ");
		
	// determine if any of the nodes are already in the node_inds
	// TODO: could use a set or other STL stuff here instead
	bool flags[4] = { false, false, false, false };
	
	for (int i=0; i < *nNode; i++) {
		if ( node_inds[i] == DT[tt].p[NODE_A] ) flags[0] = true;
		if ( node_inds[i] == DT[tt].p[NODE_B] ) flags[1] = true;
		if ( node_inds[i] == DT[tt].p[NODE_C] ) flags[2] = true;
		if ( node_inds[i] == DT[tt].p[NODE_D] ) flags[3] = true;
	}
		
	// add nodes and tetra to processing lists
	if (!flags[0]) {
		node_inds[*nNode] = DT[tt].p[NODE_A];
		*nNode += 1;
	}
	if (!flags[1]) {
		node_inds[*nNode] = DT[tt].p[NODE_B];
		*nNode += 1;
	}
	if (!flags[2]) {
		node_inds[*nNode] = DT[tt].p[NODE_C];
		*nNode += 1;
	}
	if (!flags[3]) {
		node_inds[*nNode] = DT[tt].p[NODE_D];
		*nNode += 1;
	}
	
	tet_inds[*nTet] = tt;
	*nTet += 1;
	
	IF_DEBUG(cout << "   addTet [" << tt << "] new nTet = " << *nTet << " nNode = " << *nNode 
								<< " (" << DT[tt].p[NODE_A] << " " << DT[tt].p[NODE_B] << " " << DT[tt].p[NODE_C] 
								<< " " << DT[tt].p[NODE_D] << ")" << endl);
	
	// recursively add neighbors of this tetra if needed
	if ( needTet(DT[tt].t[0], pp, tet_inds, nTet) ) 
		addTet( DT[tt].t[0], pp, node_inds, tet_inds, nNode, nTet );
	if ( needTet(DT[tt].t[1], pp, tet_inds, nTet) )
		addTet( DT[tt].t[1], pp, node_inds, tet_inds, nNode, nTet );
	if ( needTet(DT[tt].t[2], pp, tet_inds, nTet) )
		addTet( DT[tt].t[2], pp, node_inds, tet_inds, nNode, nTet );
	if ( needTet(DT[tt].t[3], pp, tet_inds, nTet) )
		addTet( DT[tt].t[3], pp, node_inds, tet_inds, nNode, nTet );
}

double ArepoMesh::ccVolume(double *ci, double *cj, double *ck, double *ct)
{
	double xt = ct[0],    yt = ct[1],    zt = ct[2];
	double xi = ci[0]-xt, yi = ci[1]-yt, zi = ci[2]-zt;
	double xj = cj[0]-xt, yj = cj[1]-yt, zj = cj[2]-zt;
	double xk = ck[0]-xt, yk = ck[1]-yt, zk = ck[2]-zt;
	
	return xi*(yj*zk-yk*zj) + yi*(zj*xk-zk*xj) + zi*(xj*yk-xk*yj);
}
#endif

// interpolate scalar fields at position pt inside Voronoi cell SphP_ID (various methods)
int ArepoMesh::subSampleCell(const Ray &ray, Vector &pt, float *vals, int taskNum)
{
#ifdef USE_DC_CONNECTIVITY
	int dpInd = -1; // not used
	int sphInd = ray.index;
#else // USE_ALTERNATIVE_CONNECTIVITY
	int dpInd = ray.index;
	int sphInd = getSphPID(DP[ray.index].index);
#endif
	
	// check degenerate point in R3, immediate return
	if (fabs(pt.x - P[sphInd].Pos[0]) <= INSIDE_EPS &&
			fabs(pt.y - P[sphInd].Pos[1]) <= INSIDE_EPS &&
			fabs(pt.z - P[sphInd].Pos[2]) <= INSIDE_EPS)
	{
			vals[TF_VAL_DENS]   = SphP[sphInd].Density;
			vals[TF_VAL_UTHERM] = SphP[sphInd].Utherm;
			return 1;
	}
				
	// zero vals we will override in this function
	vals[TF_VAL_DENS]   = 0.0;
	vals[TF_VAL_UTHERM] = 0.0;			
				
#if defined(NATURAL_NEIGHBOR_IDW) || defined(NATURAL_NEIGHBOR_SPHKERNEL)

/* exponent of distance, greater values assign more influence to values closest to the 
 * interpolating point, approaching piecewise constant for large POWER_PARAM.
 * in N dimensions, if p <= N, the interpolated values are dominated by points far away,
 * which is rather bizarre. note: p is actually 2p since we skip the sqrt.
 */
#define POWER_PARAM 2.0

/* use <1 for more smoothing, makes hsml_used bigger than hsml_ngb_max
 * use >1 for less smoothing, make hsml_used smaller than hsml_ngb_max 
 *   (at least some some neighbors will not contribute, and if this is too big, the
 *    containing cell may also not contribute, leading to black holes)
 */
#define HSML_FAC 0.95

		float dx,dy,dz,xtmp,ytmp,ztmp;
		float weight, weightsum=0, distsq;

#ifdef NATURAL_NEIGHBOR_SPHKERNEL
		// for SPHKERNEL first need to pick smoothing length
		float hinv = HSML_FAC * calcNeighborHSML(sphInd,pt,dpInd);
#endif
		
#ifndef BRUTE_FORCE

#ifdef USE_DC_CONNECTIVITY

		int edge = SphP[sphInd].first_connection;
		int last_edge = SphP[sphInd].last_connection;
		
#ifdef NATURAL_NEIGHBOR_INNER
		int pri_neighbor_inds[20];
#endif
		int k=0;
		
		edge = SphP[sphInd].first_connection;
		
		while(edge >= 0) {
			int sphp_neighbor = DC[edge].index;
			
			// could connect to bounding tetra if we don't have a ghost across this boundary
			if ( sphp_neighbor < 0 
#ifdef NO_GHOST_CONTRIBS
                        || DC[edge].dp_index >= NumGas
#endif
			) {
				if ( edge == last_edge )
					break;
				edge = DC[edge].next;
				continue;
			}
			
#ifdef NATURAL_NEIGHBOR_INNER
			pri_neighbor_inds[k++] = sphp_neighbor;
		
			// loop over neighbors of neighbors
			int inner_edge = SphP[sphp_neighbor].first_connection;
			int inner_last_edge = SphP[sphp_neighbor].last_connection;
			
			while(inner_edge >= 0) {
				// could connect to bounding tetra if we don't have a ghost across this boundary
				if ( DC[inner_edge].index < 0 ) {
					if ( inner_edge == inner_last_edge )
						break;
					inner_edge = DC[inner_edge].next;
					continue;
				}
				
				// skip any neighbors we've already done
				for( int i=0; i < k; i++ ) {
				  if ( DC[inner_edge].index == pri_neighbor_inds[i] ) {
						if(inner_edge == inner_last_edge)
							break;

						inner_edge = DC[inner_edge].next;
						continue;
					}
				}
				
				// neighbor of neighbor IDW
				dx = NGB_PERIODIC_LONG_X(P[ DC[inner_edge].index ].Pos[0] - pt.x);
				dy = NGB_PERIODIC_LONG_Y(P[ DC[inner_edge].index ].Pos[1] - pt.y);
				dz = NGB_PERIODIC_LONG_Z(P[ DC[inner_edge].index ].Pos[2] - pt.z);
				distsq = dx*dx + dy*dy + dz*dz;
			  
#ifdef NATURAL_NEIGHBOR_IDW
				weight = 1.0 / pow( sqrtf(distsq),POWER_PARAM );
#else // NATURAL_NEIGHBOR_SPHKERNEL
				weight = sph_kernel( sqrtf(distsq),hinv );
#endif
				weightsum += weight;
				
				vals[TF_VAL_DENS]   += SphP[ DC[inner_edge].index ].Density * weight;
				vals[TF_VAL_UTHERM] += SphP[ DC[inner_edge].index ].Utherm * weight;
			
			  if(inner_edge == inner_last_edge)
					break;

				inner_edge = DC[inner_edge].next;
			}
#endif // NATURAL_NEIGHBOR_INNER	
		
			dx = NGB_PERIODIC_LONG_X(P[sphp_neighbor].Pos[0] - pt.x);
			dy = NGB_PERIODIC_LONG_Y(P[sphp_neighbor].Pos[1] - pt.y);
			dz = NGB_PERIODIC_LONG_Z(P[sphp_neighbor].Pos[2] - pt.z);
			distsq = dx*dx + dy*dy + dz*dz;
			  
#ifdef NATURAL_NEIGHBOR_IDW
			weight = 1.0 / pow( sqrtf(distsq),POWER_PARAM );
#else // NATURAL_NEIGHBOR_SPHKERNEL
			weight = sph_kernel( sqrtf(distsq),hinv );
#endif
			weightsum += weight;
				
			vals[TF_VAL_DENS]   += SphP[ sphp_neighbor ].Density * weight;
			vals[TF_VAL_UTHERM] += SphP[ sphp_neighbor ].Utherm * weight;
				
			// move to next neighbor
			if(edge == last_edge)
				break;
				
#ifdef DEBUG
			if (DC[edge].next == edge || DC[edge].next < 0 || sphp_neighbor < 0)
			  terminate(" what is going on ");
#endif				

			edge = DC[edge].next;
		}
		
#else // USE_ALTERNATIVE_CONNECTIVITY

#ifdef NATURAL_NEIGHBOR_INNER
		terminate("Error: INNER for alt connectivity not implemented.");
#endif

		const int start_edge = midpoint_idx[dpInd].first;
		const int n_edges    = midpoint_idx[dpInd].second;
		
		for ( int k=0; k < n_edges; k++ ) {
			int sphp_neighbor = getSphPID( opposite_points[start_edge+k] );
			
			dx = NGB_PERIODIC_LONG_X(P[sphp_neighbor].Pos[0] - pt.x);
			dy = NGB_PERIODIC_LONG_Y(P[sphp_neighbor].Pos[1] - pt.y);
			dz = NGB_PERIODIC_LONG_Z(P[sphp_neighbor].Pos[2] - pt.z);
			distsq = dx*dx + dy*dy + dz*dz;
			  
#ifdef NATURAL_NEIGHBOR_IDW
			weight = 1.0 / pow( sqrtf(distsq),POWER_PARAM );
#else // NATURAL_NEIGHBOR_SPHKERNEL
			weight = sph_kernel( sqrtf(distsq),hinv );
#endif
			weightsum += weight;
				
			vals[TF_VAL_DENS]   += SphP[ sphp_neighbor ].Density * weight;
			vals[TF_VAL_UTHERM] += SphP[ sphp_neighbor ].Utherm * weight;
		}
		
#endif // CONNECTIVITY

#ifdef NO_GHOST_CONTRIBS
		if( sphInd < NumGas )
		{
#endif
		// add in primary parent
		dx = NGB_PERIODIC_LONG_X(P[sphInd].Pos[0] - pt.x);
		dy = NGB_PERIODIC_LONG_Y(P[sphInd].Pos[1] - pt.y);
		dz = NGB_PERIODIC_LONG_Z(P[sphInd].Pos[2] - pt.z);
		distsq = dx*dx + dy*dy + dz*dz;
			  
#ifdef NATURAL_NEIGHBOR_IDW
		weight = 1.0 / pow( sqrtf(distsq),POWER_PARAM );
#else // NATURAL_NEIGHBOR_SPHKERNEL
		weight = sph_kernel( sqrtf(distsq),hinv );
#endif
		weightsum += weight;
		
		vals[TF_VAL_DENS]   += SphP[sphInd].Density * weight;
		vals[TF_VAL_UTHERM] += SphP[sphInd].Utherm * weight;
#ifdef NO_GHOST_CONTRIBS
		}
#endif
		
#else // BRUTE_FORCE

		// brute force loop over NumGas
		for( int sphp_neighbor = 0; sphp_neighbor < NumGas; sphp_neighbor++ ) {			
			dx = NGB_PERIODIC_LONG_X(P[sphp_neighbor].Pos[0] - pt.x);
			dy = NGB_PERIODIC_LONG_Y(P[sphp_neighbor].Pos[1] - pt.y);
			dz = NGB_PERIODIC_LONG_Z(P[sphp_neighbor].Pos[2] - pt.z);
			distsq = dx*dx + dy*dy + dz*dz;
			  
#ifdef NATURAL_NEIGHBOR_IDW
			weight = 1.0 / pow( sqrtf(distsq),POWER_PARAM );
#else // NATURAL_NEIGHBOR_SPHKERNEL
			hinv = 4.0;
			weight = sph_kernel( sqrtf(distsq),hinv );
#endif
			weightsum += weight;
				
			vals[TF_VAL_DENS]   += SphP[ sphp_neighbor ].Density * weight;
			vals[TF_VAL_UTHERM] += SphP[ sphp_neighbor ].Utherm * weight;
		}
			
#endif // BRUTE_FORCE

		// normalize weights
		vals[TF_VAL_DENS]   /= weightsum;
		vals[TF_VAL_UTHERM] /= weightsum;
			
#endif // NATURAL_NEIGHBOR_IDW or NATURAL_NEIGHBOR_SPHKERNEL

/* -------------------------------------------------------------------------------------- */

#ifdef NATURAL_NEIGHBOR_INTERP
		int tlast = 0;
		float weight, weightsum = 0.0;
		
		double dp_old_vol[AUXMESH_ALLOC_SIZE/2];
		double dp_new_vol[AUXMESH_ALLOC_SIZE/2];

		// recreate the voronoi cell of the parent's neighbors (auxiliary mesh approach):
		init_clear_auxmesh(&AuxMeshes[taskNum]);
			
		// construct new auxiliary mesh around pt
#ifdef USE_DC_CONNECTIVITY
		int edge = SphP[sphInd].first_connection;
		int last_edge = SphP[sphInd].last_connection;
		
		while(edge >= 0) {
			int dp_neighbor = DC[edge].dp_index;
			
			// could connect to bounding tetra if we don't have a ghost across this boundary
			if ( DC[edge].index < 0 ) {
				if ( edge == last_edge )
					break;
				edge = DC[edge].next;
				continue;
			}
#else // USE_ALTERNATIVE_CONNECTIVITY
		const int start_edge = midpoint_idx[dpInd].first;
		const int n_edges    = midpoint_idx[dpInd].second;
		
		for( int k = 0; k < n_edges; k++ ) {
			int dp_neighbor = opposite_points[start_edge+k];
#endif
		
			if(AuxMeshes[taskNum].Ndp + 2 >= AuxMeshes[taskNum].MaxNdp)
				terminate("AuxMesh for NNI exceeds maximum size.");
				
			// insert point
			AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp] = Mesh.DP[ dp_neighbor ];
			
			// wrap this DP point to near our sample point if necessary
			// this doesn't necessary work, if some neighbors are wrapped and others aren't we will not get the correct
			// periodic volume change? in any case, only matters with ghost neighbors, deal with later
			//periodic_wrap_DP_point( AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp], pt );
				
			IF_DEBUG(cout << "   insertN (orig x=" << Mesh.DP[ dp_neighbor ].x << " y=" << Mesh.DP[ dp_neighbor ].y
			              << " z=" << Mesh.DP[ dp_neighbor ].z << ") (wrapped x=" 
										<< AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp].x 
										<< " y=" << AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp].y 
										<< " z=" << AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp].z 
										<< ") tlast=" << setw(3) << tlast << " totnum=" << setw(2) << AuxMeshes[taskNum].Ndp+1 << endl);				
				
		  //set_integers_for_point(&AuxMeshes[taskNum], AuxMeshes[taskNum].Ndp);
			set_integers_for_pointer( &AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp] );
		  tlast = insert_point_new(&AuxMeshes[taskNum], AuxMeshes[taskNum].Ndp, tlast);
		  AuxMeshes[taskNum].Ndp++;
				
#ifdef USE_DC_CONNECTIVITY
			// move to next neighbor
			if(edge == last_edge)
				break;
				
#ifdef DEBUG		 
			if (DC[edge].next == edge || DC[edge].next < 0)
			  terminate(" what is going on ");
#endif
				
			edge = DC[edge].next;
#endif
		}
		
		// add primary parent as well
#ifdef USE_DC_CONNECTIVITY
		AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp] = Mesh.DP[sphInd];
#else
		AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp] = Mesh.DP[dpInd];
#endif

		// wrap
		//periodic_wrap_DP_point( AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp], pt);
		
		IF_DEBUG(cout << "   insertP (orig x=" << Mesh.DP[ sphInd ].x << " y=" << Mesh.DP[ sphInd ].y
		              << " z=" << Mesh.DP[ sphInd ].z << ") (wrapped x=" 
									<< AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp].x 
									<< " y=" << AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp].y 
									<< " z=" << AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp].z 
									<< ") tlast=" << setw(3) << tlast << " totnum=" << setw(2) << AuxMeshes[taskNum].Ndp+1 << endl);		
		
		set_integers_for_pointer( &AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp] );
		tlast = insert_point_new(&AuxMeshes[taskNum], AuxMeshes[taskNum].Ndp, tlast);
		AuxMeshes[taskNum].Ndp++;
		
		// compute old circumcircles and volumes
		compute_circumcircles(&AuxMeshes[taskNum]);
		compute_auxmesh_volumes(&AuxMeshes[taskNum], dp_old_vol);
		
#ifdef DEBUG
		cout << "   old volumes:";
		for (int k = 0; k < AuxMeshes[taskNum].Ndp; k++) {
				cout << " [" << k << "] " << dp_old_vol[k];
		}
		cout << endl;
#endif		
		
		// add interpolate point
		AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp].x = pt.x;
		AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp].y = pt.y;
		AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp].z = pt.z;
		
		set_integers_for_pointer( &AuxMeshes[taskNum].DP[AuxMeshes[taskNum].Ndp] );
		tlast = insert_point_new(&AuxMeshes[taskNum], AuxMeshes[taskNum].Ndp, tlast);
		AuxMeshes[taskNum].Ndp++;
		
		IF_DEBUG(cout << "   inserted interpolate new tlast=" << tlast << " totnum=" << AuxMeshes[taskNum].Ndp << endl);
	
		// compute new circumcircles and volumes
		compute_circumcircles(&AuxMeshes[taskNum]);
		compute_auxmesh_volumes(&AuxMeshes[taskNum], dp_new_vol);

#ifdef DEBUG
		cout << "   new volumes:";
		for (int k = 0; k < AuxMeshes[taskNum].Ndp; k++) {
				cout << " [" << k << "] " << dp_new_vol[k];
		}
		cout << endl;
#endif

		// calculate scalar value based on neighbor values and area fraction weights
#ifdef USE_DC_CONNECTIVITY
		edge = SphP[sphInd].first_connection;
		int k = 0;
		
		while(edge >= 0) {
			int sphp_neighbor = DC[edge].index;
			// could connect to bounding tetra if we don't have a ghost across this boundary
			if ( sphp_neighbor < 0 ) {
				if ( edge == last_edge )
					break;
				edge = DC[edge].next;
				continue;
			}
#else // USE_ALTERNATIVE_CONNECTIVITY
		for ( int k=0; k < n_edges; k++ ) {
			int sphp_neighbor = getSphPID( DP[opposite_points[start_edge+k]].index );
#endif

			weight = dp_old_vol[k] - dp_new_vol[k];
			weightsum += weight;
			
			vals[TF_VAL_DENS]   += SphP[sphp_neighbor].Density * weight;
			vals[TF_VAL_UTHERM] += SphP[sphp_neighbor].Utherm * weight;
			
#ifdef USE_DC_CONNECTIVITY
			k++;
			
			// move to next neighbor
			if(edge == last_edge)
				break;
				
#ifdef DEBUG		 
			if (DC[edge].next == edge || DC[edge].next < 0)
			  terminate(" what is going on ");
#endif

			edge = DC[edge].next;	
#endif			
		}
		
		// add primary parent
		weight = dp_old_vol[AuxMeshes[taskNum].Ndp-2] - dp_new_vol[AuxMeshes[taskNum].Ndp-2];
		weightsum += weight;
		
		vals[TF_VAL_DENS]   += SphP[sphInd].Density * weight;
		vals[TF_VAL_UTHERM] += SphP[sphInd].Utherm * weight;

		IF_DEBUG(cout << "   weightsum = " << weightsum << " inserted vol = " << dp_new_vol[AuxMeshes[taskNum].Ndp-1] << endl);
		
		// do the volumes lost by all the natural neighbors add up to the sample pt cell volume?
		//if ( fabs(weightsum - dp_new_vol[AuxMeshes[taskNum].Ndp-1]) / dp_new_vol[AuxMeshes[taskNum].Ndp-1] > 1e-1 )
		//  terminate("NNI weight error is large.");
		
		// normalize by volume of last added cell (around interp point)
		vals[TF_VAL_DENS]   /= dp_new_vol[AuxMeshes[taskNum].Ndp-1];
		vals[TF_VAL_UTHERM] /= dp_new_vol[AuxMeshes[taskNum].Ndp-1];

#endif // NATURAL_NEIGHBOR_INTERP

/* -------------------------------------------------------------------------------------- */

#ifdef DTFE_INTERP
		// use DT[tt].p[0] as the x_i sample point
		int tt0_DPID   = DT[ray.tetra].p[0];
		
		// if p[0] is part of the bounding tetra (we are on the edge), cannot do DTFE
		if (tt0_DPID < 0) {
			IF_DEBUG(cout << "  dtfe: p[0] = " << tt0_DPID << ", skipping!" << endl);
		  return 0;
		}
		
		int tt0_SphPID = getSphPID(DP[tt0_DPID].index);
		
		Vector p0cen( DP[tt0_DPID].x, DP[tt0_DPID].y, DP[tt0_DPID].z );
		Vector tetraGrad( DT_grad[3*ray.tetra+0], DT_grad[3*ray.tetra+1], DT_grad[3*ray.tetra+2] );
		
		pt -= p0cen; // make relative to p[0] position (not Voronoi center)

		// apply the (linear) gradient to the sampling point
		vals[TF_VAL_DENS]   = SphP[tt0_SphPID].Density + Dot(tetraGrad,pt);
		vals[TF_VAL_UTHERM] = SphP[sphInd].Utherm; // TODO
		
#ifdef DEBUG
		cout << "  dtfe: tt = " << ray.tetra << " p0 = " << tt0_DPID << " sph0 = " << tt0_SphPID << endl;
		pt.print("  dtfe: relative point ");
		cout << "  dtfe: dens = " << vals[TF_VAL_DENS] << " sphDens = " << SphP[tt0_SphPID].Density << endl;
#endif
		
#endif // DTFE

/* -------------------------------------------------------------------------------------- */

#ifdef NNI_WATSON_SAMBRIDGE
	// const int access_tri_ws[4][3] = { {1, 2, 3}, {2, 3, 0},	{3, 0, 1}, {0, 1, 2} };
	// const int access_nodes_ws[4] = {3,0,1,2}; // used via NODE_X
	
	// circumcenters of five new tetras composed of 3 original vertices and the pt
	double cca[3]; // pbcd
	double ccb[3]; // padc
	double ccc[3]; // pabd
	double ccd[3]; // pacb
	double cct[3]; // abcd, the original tetra
	
	// volumes of four new tetras each composed of 1 original tetra face and the pt
	double va,vb,vc,vd;
	
	// total volume sum of voronoi cell that would be added (pieces are in DP_vols)
	double vol_sum = 0.0;
	
	// double triplet for our interp point
	point pp;
	pp.x = pt.x;
	pp.y = pt.y;
	pp.z = pt.z;
	set_integers_for_pointer(&pp);
	
	// compute the tetra/node lists corresponding to natural neighbors
	int nNode = 0;
	int nTet = 0;
	
	int tet_inds[MAX_NUM_TETS];
	int node_inds[MAX_NUM_NODES];
	
	// start the recursive search with the current parent tetra
	ArepoMesh::addTet(ray.tetra, &pp, &node_inds[0], &tet_inds[0], &nNode, &nTet);
	
	// zero volumes we will accumulate into and make nodeIndex -> DP_vols_index map
	float DP_vols[MAX_NUM_NODES];
	multimap<int,int> mm;
	
	for ( int i=0; i < nNode; i++ ) {
		DP_vols[ i ] = 0.0;
		mm.insert( std::pair<int,int>( node_inds[i], i ) );
	}
		
	// loop over natural neighbor tetras
	bool okFlag = true;
	
	for ( int i=0; i < nTet; i++ )
	{
		int tet = tet_inds[i];
		
		// compute five new circumcenters
		// cyclic ordering correct? (0=a, 1=b, 2=c, 3=d is incorrect)
		// (Arepo orientation: 3 lies above the oriented triangle formed by 0,1,2)
		// --> (3=a, 0=b, 1=c, 2=d if tri is CCW) (3=a, 0=b, 2=c, 1=d if tri is CW)
		okFlag &= calc_circumcenter(T, &pp, DT[tet].p[NODE_B], DT[tet].p[NODE_C], DT[tet].p[NODE_D], &cca[0]);
		okFlag &= calc_circumcenter(T, &pp, DT[tet].p[NODE_A], DT[tet].p[NODE_D], DT[tet].p[NODE_C], &ccb[0]);
		okFlag &= calc_circumcenter(T, &pp, DT[tet].p[NODE_A], DT[tet].p[NODE_B], DT[tet].p[NODE_D], &ccc[0]);
		okFlag &= calc_circumcenter(T, &pp, DT[tet].p[NODE_A], DT[tet].p[NODE_C], DT[tet].p[NODE_B], &ccd[0]);
		okFlag &= calc_circumcenter(T, &(DP[DT[tet].p[NODE_A]]), 
		                          DT[tet].p[NODE_B], DT[tet].p[NODE_C], DT[tet].p[NODE_D], &cct[0]);
		
		// compute four volumes
		va = ArepoMesh::ccVolume(&ccb[0],&ccc[0],&ccd[0],&cct[0]);
		vb = ArepoMesh::ccVolume(&cca[0],&ccd[0],&ccc[0],&cct[0]);
		vc = ArepoMesh::ccVolume(&cca[0],&ccb[0],&ccd[0],&cct[0]);
		vd = ArepoMesh::ccVolume(&cca[0],&ccc[0],&ccb[0],&cct[0]);
		
		// accumulate volumes onto vertices (nodes)
		DP_vols[ mm.find( DT[tet].p[NODE_A] )->second ] += va;
		DP_vols[ mm.find( DT[tet].p[NODE_B] )->second ] += vb;
		DP_vols[ mm.find( DT[tet].p[NODE_C] )->second ] += vc;
		DP_vols[ mm.find( DT[tet].p[NODE_D] )->second ] += vd;
		
		vol_sum += va + vb + vc + vd;
		
		IF_DEBUG(cout << "   tet [" << setw(2) << i << "] [" << setw(4) << tet << "] va = " 
									<< va << " vb = " << vb << " vc = " 
									<< vc << " vd = " << vd << " vol_sum = " << vol_sum << endl);
	}
	
	IF_DEBUG(cout << "   final vol_sum = " << vol_sum << endl);
	
	if ( !okFlag ) {
		IF_DEBUG(cout << "  errFlag = 1, poor circumcenters, skipping point!" << endl);
		return 0;
	}
	
	// we have the final volumes and the normalization, loop once more over the contributing tetras
	double wt_total = 0.0; // debugging

	for ( int i=0; i < nNode; i++ )
	{
	  int sphInd = getSphPID(DP[node_inds[i]].index);
		vals[TF_VAL_DENS]   += SphP[sphInd].Density * DP_vols[ node_inds[i] ];
		vals[TF_VAL_UTHERM] += SphP[sphInd].Utherm  * DP_vols[ node_inds[i] ];
		
		IF_DEBUG(cout << "   add node [i " << setw(2) << i << "] [sphInd " << setw(3) << sphInd
									<< "] Dens = " << SphP[sphInd].Density 
									<< " wt = " << DP_vols[ mm.find( node_inds[i] )->second ]/vol_sum 
									<< " ( " << SphP[sphInd].Density * DP_vols[ mm.find( node_inds[i] )->second ] << " )" << endl);
									
		wt_total += DP_vols[ mm.find( node_inds[i] )->second ];
	}
	
	IF_DEBUG(cout << "   wt_total = " << wt_total << " ( error = " << wt_total - vol_sum << " )" << endl);
	
	if ( vol_sum <= 0.0 )
	  terminate("Watson-Sambridge produced negative volume for inserted Vcell.");
	//if ( fabs(wt_total - vol_sum)/vol_sum > 1e-1 )
	//  terminate("Watson-Sambridge weights failed to sum up to total volume (%g %g %g).",pt.x,pt.y,pt.z);
	
	vals[TF_VAL_DENS]   /= vol_sum;
	vals[TF_VAL_UTHERM] /= vol_sum;

#endif

/* -------------------------------------------------------------------------------------- */

#ifdef NNI_LIANG_HALE
		terminate("TODO");
#endif

/* -------------------------------------------------------------------------------------- */

#ifdef CELL_GRADIENTS_ONLY
		// old:
    pt -= Vector( P[sphInd].Pos );	// make relative to cell center	
		
		// Voronoi stencil-based linear gradient
		Vector sphDensGrad( SphP[sphInd].Grad.drho );
		
		vals[TF_VAL_DENS]   = SphP[sphInd].Density + Dot(sphDensGrad,pt);
		vals[TF_VAL_UTHERM] = SphP[sphInd].Utherm; // no utherm gradient unless MATERIALS
		
#endif // CELL_GRADIENTS_ONLY

		return 1;
}

#endif // ENABLE_AREPO
