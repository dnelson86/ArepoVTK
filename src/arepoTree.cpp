/*
 * arepoTree.cpp
 * dnelson
 */
 
#include "transform.h"
#include "spectrum.h"
#include "volume.h"
#include "transfer.h"
#include "arepo.h"
#include "util.h" // for numberOfCores()

#define NGB_LIST_FAC 4 // times requested nNGB safety margin

ArepoTree::ArepoTree(const TransferFunction *tf)
{
		IF_DEBUG(cout << "ArepoTree() constructor." << endl);
		
		// transfer function and sampling setup
		transferFunction = tf;
		viStepSize       = Config.viStepSize;
		
		sampleWt = 1.0f; //All.BoxSize / pow(NumGas,0.333);
		
		if (viStepSize)
			sampleWt *= viStepSize;
		
		// set pointers into Arepo data structures
		
		if (Config.verbose)
				cout << "[" << ThisTask << "] ArepoMesh: NumGas = " << NumGas << " NumPart = " << NumPart << endl << endl;
		
		// allocate variable size tree search return list
		//int numLists = numberOfCores();
		//for(int i = 0; i < numLists; i++)
		//	varNGBLists.push_back( new int[Config.nTreeNGB * NGB_LIST_FAC] );
			
		//varNGBList = new int[Config.nTreeNGB * NGB_LIST_FAC]; //OLD
		
		// boxsize
		extent = BBox(Point(0.0,0.0,0.0),Point(All.BoxSize,All.BoxSize,All.BoxSize));
		
		if( !All.BoxSize ) {
			cout << "Error: All.BoxSize=0, likely structure mismatch." << endl;
			exit(1180);
		}		
		
		IF_DEBUG(extent.print(" ArepoTree extent "));		
}

ArepoTree::~ArepoTree()
{
		// de-allocate variable neighbor lists
		//int numLists = numberOfCores();
		//for(int i = 0; i < numLists; i++)
		//	delete varNGBLists[i];
}

bool ArepoTree::FindNeighborList(Point &pt, float hsml, int *numngb_int, vector<float> &vals)
{
	int numngb = 0;
#if defined(NATURAL_NEIGHBOR_SPHKERNEL) || defined(NATURAL_NEIGHBOR_IDW)
	// based on ngb_treefind_variable (e.g. CalcTHValWt)
	int node, p;
	struct NgbNODE *current;
	double dx, dy, dz, xtmp, ytmp, ztmp, r, wk, weight, h2, r2;
	float search_min[3], search_max[3], search_max_Lsub[3], search_min_Ladd[3];
	
#ifdef NATURAL_NEIGHBOR_SPHKERNEL
	double u, hinv, hinv3;
	
	hinv = 1.0 / hsml;
	hinv3 = hinv * hinv * hinv;
#endif
	
	h2 = hsml * hsml;
	r2 = 0.0;
	weight = 0.0;
	
	// starting node
	node = Ngb_MaxPart;
		
	// search bounds
	search_min[0] = pt.x - hsml;
	search_min[1] = pt.y - hsml;
	search_min[2] = pt.z - hsml;
	search_max[0] = pt.x + hsml;
	search_max[1] = pt.y + hsml;
	search_max[2] = pt.z + hsml;

	search_max_Lsub[0] = search_max[0] - boxSize_X;
	search_max_Lsub[1] = search_max[1] - boxSize_Y;
	search_max_Lsub[2] = search_max[2] - boxSize_Z;

	search_min_Ladd[0] = search_min[0] + boxSize_X;
	search_min_Ladd[1] = search_min[1] + boxSize_Y;
	search_min_Ladd[2] = search_min[2] + boxSize_Z;
	
	// start search
	while(node >= 0)
	{
		if(node < Ngb_MaxPart)  // single particle
		{
			p = node;
			node = Ngb_Nextnode[node];

			if(P[p].Type > 0) // not gas particle
				continue;
			
			dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - pt.x);
			if(dx > hsml)
				continue;
				
			dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - pt.y);
			if(dy > hsml)
				continue;
				
			dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - pt.z);
			if(dz > hsml)
				continue;
	    
			r2 = dx * dx + dy * dy + dz * dz;
			
			if( r2 < h2 )
			{
				//if( numngb >= Config.nTreeNGB * NGB_LIST_FAC )
				//	terminate("going to buffer overrun, increase NGB_LIST_FAC");
					
				r = sqrt(r2);
#ifdef NATURAL_NEIGHBOR_IDW
				wk = 1.0 / pow( r,POWER_PARAM );
#else

#ifdef NATURAL_NEIGHBOR_SPHKERNEL
				u = r * hinv;
				
				if(u < 0.5)
					wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
				else
					wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

				wk = NORM_COEFF * wk / hinv3;
#endif

#endif // NATURAL_NEIGHBOR_IDW

				// accumulate weights and add to neighbor list
				weight += wk;
				
				//varNGBList[numngb++] = p;
				numngb++;
				
				// TEST: store directly weighted values
				addValsContribution( vals, p, wk );
			}
			
		}
		else if(node < Ngb_MaxPart + Ngb_MaxNodes) // internal node
		{
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
			node = Ngb_Nextnode[node - Ngb_MaxNodes];
			continue;
		}
	}

	// normalize vals
	weight = 1.0 / weight;
	
	for( unsigned int i=0; i < vals.size(); i++ )
		vals[i] *= weight;
	
#ifdef DEBUG
	cout << "FindNeighborList(): numngb = " << numngb << " weight = " << weight 
	     << " ( dens = " << vals[TF_VAL_DENS] << " utherm = " << vals[TF_VAL_TEMP] << " ) " << endl;
#endif	
#endif // SPH/IDW
	*numngb_int = numngb;
	
	if( !numngb )
	  return false; // skip

	return true;
}

bool ArepoTree::AdvanceRayOneStep(const Ray &ray, double *t0, double *t1, 
																	Spectrum &Lv, Spectrum &Tr, int threadNum)
{
		double min_t_old, min_t_new;
		vector<float> vals(TF_NUM_VALS, 0.0);
		int numngb_int;
		bool status;
		
		// verify task
		if (ray.task != ThisTask)
			terminate("Ray on wrong task.");
	
		// setup stepping: strict in world space
		min_t_old = *t0 + ray.depth * viStepSize;
		min_t_new = *t0 + (ray.depth+1) * viStepSize;
		Point prev_sample_pt( ray(min_t_old) );
	
		// check if exiting box and failed to exit a face
		if ( !extent.Inside( prev_sample_pt ) ) {
			// set intersection with box face to allow for final contribution to ray
			IF_DEBUG(cout << " current ray pos outside box, ok!" << endl);
			min_t_new = ray.max_t;
		}
		
		IF_DEBUG(prev_sample_pt.print("  prev_sample_pt "));
					
		// sample point
		min_t_new = Clamp(min_t_new,*t0,*t1); // clamp min_t_new to avoid integrating outside the box
		Point midpt( ray(min_t_new) );
			
		// tree search for N nearest neighbors
		status = FindNeighborList( midpt, ray.prevHSML, &numngb_int, vals );
		
		// adjust hsml (for the next sample point) based on difference between requested nTreeNGB
		// and the number of neighbors found with this current hsml
		float newHsml = (Config.nTreeNGB - numngb_int);
		newHsml *= 0.2 / Config.nTreeNGB;
		newHsml += 1.0f;
			
		newHsml = Clamp(newHsml, 0.5, 2.0); // change by at most a factor of two per sample point
		ray.prevHSML = newHsml * ray.prevHSML;
#ifdef DEBUG
		cout << " newHsml = " << setprecision(6) << newHsml << " prevHSML = " << setprecision(6) << ray.prevHSML << endl;
#endif

		// sample quantities at this position (now have interpolated densisty,temp,etc)
		// i.e. fill vals vector
		// currently done in FindNeighborList()
		
		if(status)
		{
			// optical depth
			Spectrum stepTau(0.0);
			stepTau += transferFunction->sigma_t() * ( vals[TF_VAL_DENS] ) * viStepSize;
			//Tr *= Exp(-stepTau); // reduce transmittance for optical depth
						
			// compute emission-only source term using transfer function
			if(status)
				Lv += Tr * transferFunction->Lve(vals) * sampleWt;
		}

		// update ray: transfer to next voronoi cell (possibly on different task)
		ray.depth++;
		ray.task = 0;
		ray.min_t = Clamp(min_t_new,ray.min_t,ray.max_t);
			
		IF_DEBUG(cout << " updated ray new task = " << ray.task << " depth = " << ray.depth 
									<< " min_t = " << ray.min_t << endl);
			
		if (fabs(ray.min_t - ray.max_t) <= INSIDE_EPS) {
				// apparently this ray is done?
				IF_DEBUG(cout << " min_t == t1 = " << *t1 << ", ray done." << endl);
				return false;
		}
		
		return true;		
}
