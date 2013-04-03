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

#ifdef NNI_WATSON_SAMBRIDGE
bool ArepoMesh::needTet(int tt, point *pp)
{
	// check if pp is inside tetra tt, use non-exact if possible
	int ret;
	
	ret = InSphere_Errorbound( &DP[DT[tt].p[0]], &DP[DT[tt].p[1]], &DP[DT[tt].p[2]],
														 &DP[DT[tt].p[3]], pp );
														 
	if ( ret == 0 )
		terminate(" add exact ");
		
	return ret;
}

void ArepoMesh::addTet(int tt, point *pp, int *node_inds, int *tet_inds, int *nNode, int *nTet)
{
	// add nodes and tetra to processing lists
	node_inds[*nNode+0] = DT[tt].p[0];
	node_inds[*nNode+1] = DT[tt].p[1];
	node_inds[*nNode+2] = DT[tt].p[2];
	node_inds[*nNode+3] = DT[tt].p[3];
	*nNode += 4;
	
	tet_inds[*nTet] = tt;
	*nTet += 1;
	
#ifdef DEBUG
	cout << " addTet [" << tt << "] new nTet = " << nTet << " nNode = " << nNode << " (" << DT[tt].p[0] << " " << DT[tt].p[1] << " " << DT[tt].p[2] << " " << DT[tt].p[3] << ")" << endl;
#endif
	
	// recursively add neighbors of this tetra if needed
	if ( needTet(DT[tt].t[0], pp) ) addTet( DT[tt].t[0], pp, node_inds, tet_inds, nNode, nTet );
	if ( needTet(DT[tt].t[1], pp) ) addTet( DT[tt].t[1], pp, node_inds, tet_inds, nNode, nTet );
	if ( needTet(DT[tt].t[2], pp) ) addTet( DT[tt].t[2], pp, node_inds, tet_inds, nNode, nTet );
	if ( needTet(DT[tt].t[3], pp) ) addTet( DT[tt].t[3], pp, node_inds, tet_inds, nNode, nTet );
}

double ArepoMesh::ccVolume(double *ci, double *cj, double *ck, double *ct)
{
	double xt = ct[0], yt = ct[1], zt = ct[2];
	double xi = ci[0]-xt, yi = ci[1]-yt, zi = ci[2]-zt;
	double xj = cj[0]-xt, yj = cj[1]-yt, zj = cj[2]-zt;
	double xk = ck[0]-xt, yk = ck[1]-yt, zk = ck[2]-zt;
	
	return xi*(yj*zk-yk*zj) + yi*(zj*xk-zk*xj) + zi*(xj*yk-xk*yj);
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
		float dx,dy,dz,xtmp,ytmp,ztmp;	
		float weight,weightsum=0,distsq;
		
		// loop over each neighbor
		for (int k=0; k < n_edges; k++) {
			dp_neighbor = opposite_points[start_edge + k];
			sphp_neighbor = getSphPID(DP[dp_neighbor].index);
			
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
		int edge, last_edge;
		
		// add parent to list
		int dp_neighbor, sphp_neighbor;
		float dx,dy,dz,xtmp,ytmp,ztmp;	
		float weight,weightsum=0,distsq;
		float hsml2 = 0.0;
		
		// 1. smoothing length parameter: calculate over neighbor distances
		edge = SphP[SphP_ID].first_connection;
		last_edge = SphP[SphP_ID].last_connection;
		
		while(edge >= 0) {
		  dp_neighbor = DC[edge].dp_index;
			
			dx = NGB_PERIODIC_LONG_X(DP[dp_neighbor].x - pt.x);
			dy = NGB_PERIODIC_LONG_Y(DP[dp_neighbor].y - pt.y);
			dz = NGB_PERIODIC_LONG_Z(DP[dp_neighbor].z - pt.z);
			
			distsq = dx*dx + dy*dy + dz*dz;
                                           
      if(distsq > hsml2)
        hsml2 = distsq;
				
			// move to next neighbor
			if(edge == last_edge)
				break;
				
			if (DC[edge].next == edge || DC[edge].next < 0 || dp_neighbor < 0)
			  terminate(" what is going on ");
				
			edge = DC[edge].next;
		}
		
    float hinv = HSML_FAC / sqrtf(hsml2);
         
		// 2. loop over each neighbor and add contribution
		edge = SphP[SphP_ID].first_connection;
		
		while(edge >= 0) {
		  dp_neighbor = DC[edge].dp_index;
			sphp_neighbor = getSphPID(DP[dp_neighbor].index);
			
			// calculate weight as sphkernel(distsq/h) (cubic spline kernel)
			dx = NGB_PERIODIC_LONG_X(DP[dp_neighbor].x - pt.x);
			dy = NGB_PERIODIC_LONG_Y(DP[dp_neighbor].y - pt.y);
			dz = NGB_PERIODIC_LONG_Z(DP[dp_neighbor].z - pt.z);
			
			distsq = dx*dx + dy*dy + dz*dz;
			
			weight = sph_kernel(sqrtf(distsq),hinv);
			weightsum += weight; // * SphP[SphP_ID].Volume;
			
			vals[TF_VAL_DENS]   += SphP[sphp_neighbor].Density * weight;
			vals[TF_VAL_UTHERM] += SphP[sphp_neighbor].Utherm * weight;
				
			// move to next neighbor
			if(edge == last_edge)
				break;
				
			if (DC[edge].next == edge || DC[edge].next < 0 || dp_neighbor < 0 || sphp_neighbor < 0)
			  terminate(" what is going on ");
				
			edge = DC[edge].next;
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
			sphp_neighbor = getSphPID(DP[dp_neighbor].index);
			
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
		
		int tt0_SphPID = getSphPID(DP[tt0_DPID].index);
		
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

#ifdef NNI_WATSON_SAMBRIDGE
	/* const int access_tri_ws[4][3] = {
		{1, 2, 3},
		{2, 3, 0},
		{3, 0, 1},
		{0, 1, 2}
	}; */
	
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
	pp.x = pt.y;
	pp.x = pt.z;
	set_integers_for_pointer(&pp);
	
	// compute the tetra/node lists corresponding to natural neighbors
	int nNode = 0;
	int nTet = 0;
	
	int MAX_NUM_TETS = 100;
	int MAX_NUM_NODES = 100;
	
	int tet_inds[MAX_NUM_TETS];
	int node_inds[MAX_NUM_NODES];
	
	// start the recursive search with the current parent tetra
	ArepoMesh::addTet(ray.tetra, &pp, &node_inds[0], &tet_inds[0], &nNode, &nTet);
		
	// zero volumes we will accumulate into
	for ( int i=0; i < nNode; i++ ) {
		DP_vols[node_inds[i]] = 0.0;
	}
		
	// loop over natural neighbor tetras
	for ( int i=0; i < nTet; i++ )
	{
		int tet = tet_inds[i];
		
		// compute five new circumcenters
		// cyclic ordering correct? (0=a, 1=b, 2=c, 3=d)
		calc_circumcenter(T, &pp, DT[tet].p[1], DT[tet].p[2], DT[tet].p[3], &cca[0]);
		calc_circumcenter(T, &pp, DT[tet].p[0], DT[tet].p[3], DT[tet].p[2], &ccb[0]);
		calc_circumcenter(T, &pp, DT[tet].p[0], DT[tet].p[1], DT[tet].p[3], &ccc[0]);
		calc_circumcenter(T, &pp, DT[tet].p[0], DT[tet].p[2], DT[tet].p[1], &ccd[0]);
		calc_circumcenter(T, &(DP[DT[tet].p[0]]), DT[tet].p[1], DT[tet].p[2], DT[tet].p[3], &cct[0]);
		
		// compute four volumes
		va = ArepoMesh::ccVolume(&ccb[0],&ccc[0],&ccd[0],&cct[0]);
		vb = ArepoMesh::ccVolume(&cca[0],&ccd[0],&ccc[0],&cct[0]);
		vc = ArepoMesh::ccVolume(&cca[0],&ccb[0],&ccd[0],&cct[0]);
		vd = ArepoMesh::ccVolume(&cca[0],&ccc[0],&ccb[0],&cct[0]);
		
		// accumulate volumes onto vertices (nodes)
		DP_vols[ DT[tet].p[0] ] += va;
		DP_vols[ DT[tet].p[1] ] += vb;
		DP_vols[ DT[tet].p[2] ] += vc;
		DP_vols[ DT[tet].p[3] ] += vd;
		
		vol_sum += va + vb + vc + vd;
	}
	
	// we have the final volumes and the normalization, loop once more over the contributing tetras
	for ( int i=0; i < nNode; i++ )
	{
	  int SphPID = getSphPID(DP[node_inds[i]].index);
		vals[TF_VAL_DENS]   += SphP[SphPID].Density * DP_vols[ node_inds[i] ];
		vals[TF_VAL_UTHERM] += SphP[SphPID].Utherm  * DP_vols[ node_inds[i] ];
	}
	
	vals[TF_VAL_DENS]   /= vol_sum;
	vals[TF_VAL_UTHERM] /= vol_sum;

#endif

#ifdef NNI_LIANG_HALE
		terminate("TODO");
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

#endif // ENABLE_AREPO