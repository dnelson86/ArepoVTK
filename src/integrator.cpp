/*
 * integrator.cpp
 * dnelson
 */
 
#include "integrator.h"
#include "ArepoRT.h"
#include "spectrum.h"
#include "transform.h"
#include "sampler.h"
#include "volume.h"
#include "camera.h"
#include "renderer.h"

// ------------------------------- EmissionIntegrator -------------------------------
void EmissionIntegrator::RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene)
{
    tauSampleOffset = sample->Add1D(1);
    scatterSampleOffset = sample->Add1D(1);
}

Spectrum EmissionIntegrator::Transmittance(const Scene *scene, const Renderer *renderer, const Ray &ray,
																					 const Sample *sample, RNG &rng) const
{
    if (!scene->volumeRegion) return Spectrum(1.0f);
    float step, offset;
    if (sample) {
        step = stepSize;
        offset = sample->oneD[tauSampleOffset][0];
    }
    else {
        step = 4.0f * stepSize;
        offset = rng.RandomFloat();
    }
    Spectrum tau = scene->volumeRegion->tau(ray, step, offset);
    return Exp(-tau);
}

Spectrum EmissionIntegrator::Li(const Scene *scene, const Renderer *renderer, const Ray &ray,
															  const Sample *sample, RNG &rng, Spectrum *T, int *prevEntryCell, int *prevEntryTetra, int taskNum) const
{
		IF_DEBUG(cout << "EmissionIntegrator::Li()" << endl);
		
    double t0, t1;
    if (!scene->volumeRegion || !scene->volumeRegion->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.0f) {
				IF_DEBUG(cout << " Returning! IntersectP t0 = " << t0 << " t1 = " << t1 << endl);
        *T = Spectrum(1.0f);
        return 0.0f;
    }
    // Do emission-only volume integration in _vr_
    Spectrum Lv(0.0);

    // Prepare for volume integration stepping
    int nSamples = (int)ceilf((t1-t0) / stepSize);
    float step = (t1 - t0) / nSamples;
    Spectrum Tr(1.0f);
    Point p = ray(t0), pPrev;
    Vector w = -ray.d;
    t0 += sample->oneD[scatterSampleOffset][0] * step;
		
		IF_DEBUG(cout << " t0 = " << t0 << " t1 = " << t1 << " nSamples = " << nSamples << " step = " << step << endl);
		
    for (int i = 0; i < nSamples; ++i, t0 += step) {
        // Advance to sample at _t0_ and update _T_
        pPrev = p;
        p = ray(t0);
        Ray tauRay(pPrev, p - pPrev, 0.0f, 1.0f, ray.time, ray.depth);
        Spectrum stepTau = scene->volumeRegion->tau(tauRay,0.5f * stepSize, rng.RandomFloat());
        Tr *= Exp(-stepTau);

        // Possibly terminate ray marching if transmittance is small
        if (Tr.y() < 1e-3) {
            const float continueProb = 0.5f;
            if (rng.RandomFloat() > continueProb) break;
            Tr /= continueProb;
        }

        // Compute emission-only source term at _p_
        Lv += Tr * scene->volumeRegion->Lve(p, w, ray.time);
    }
    *T = Tr;
    return Lv * step;
}


EmissionIntegrator *CreateEmissionVolumeIntegrator(const float stepSize)
{
    return new EmissionIntegrator(stepSize);
}

// ------------------------------- VoronoiIntegrator -------------------------------
void VoronoiIntegrator::Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer)
{
		// find entry voronoi cells for rays
		
		// distribute rays to appropriate start tasks
}

void VoronoiIntegrator::RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene)
{
    tauSampleOffset = sample->Add1D(1);
    scatterSampleOffset = sample->Add1D(1);
}

Spectrum VoronoiIntegrator::Transmittance(const Scene *scene, const Renderer *renderer, const Ray &ray,
																					 const Sample *sample, RNG &rng) const
{
		cout << "CHECK" << endl;
		endrun(1106);
    if (!scene->arepoMesh) return Spectrum(1.0f);
		return Spectrum(1.0f); // added just to suppress return warning, CHECK
}

Spectrum VoronoiIntegrator::Li(const Scene *scene, const Renderer *renderer, const Ray &ray,
															  const Sample *sample, RNG &rng, Spectrum *T, int *prevEntryCell, 
																int *prevEntryTetra, int threadNum) const
{
		IF_DEBUG(cout << "VoronoiIntegrator::Li()" << endl);
		
    double t0, t1;
		
    if (!scene->arepoMesh || !scene->arepoMesh->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.0f) {
				IF_DEBUG(cout << " Returning! IntersectP t0 = " << t0 << " t1 = " << t1 << endl);
        *T = Spectrum(1.0f);
        return 0.0f;
    }
			
		// do emission only volume integration in AM
		Spectrum Lv(0.0);
		Spectrum Tr(1.0f);
		
		// propagate ray to arepo box, set exit point
		ray.min_t = t0;
		ray.max_t = t1;
		
#ifdef DEBUG
		ray(t0).print(" hitbox ");
		ray(t1).print(" exitbox ");
#endif		
		
		// maximum ray integration length from config
		if (Config.rayMaxT && Config.rayMaxT < ray.max_t)
				ray.max_t = Config.rayMaxT;
				
		IF_DEBUG(cout << " t0 = " << t0 << " t1 = " << t1
									<< " ray.min_t = " << ray.min_t << " ray.max_t = " << ray.max_t << endl);		
		
		// no actual gas cells (no snapshot loaded, doing something else)?
		if( !NumGas )
			return Lv;

		// find the voronoi cell the ray will enter (or be in) first
		scene->arepoMesh->LocateEntryCell(ray, prevEntryCell);
		
#if defined(DEBUG_VERIFY_ENTRY_CELLS)
		Point pos = ray(ray.min_t);
				
		scene->arepoMesh->VerifyPointInCell(ray.index,pos);
#endif
		
#if defined(DTFE_INTERP) || defined(NNI_WATSON_SAMBRIDGE) || defined(NNI_LIANG_HALE)
		// find the delaunay tetra the ray will enter (or be in) first
		scene->arepoMesh->LocateEntryTetra(ray, prevEntryTetra);
#endif

		// TODO: exchange?
		
		// advance ray through voronoi cells
		int count = 0;
#ifdef DEBUG
		Point p = ray(ray.min_t);
		cout << " VoronoiIntegrator::Li(iter=" << count << ") Lv.y = " << setw(6) << Lv.y()
				 << " Tr.y = " << Tr.y() << " ray.x = " << setw(5) << p.x 
				 << " ray.y = " << setw(5) << p.y << " ray.z = " << setw(5) << p.z << endl;
#endif			 

		while( true ) {
				count++;
				
				// TODO: SIGTERM
				// if( sigTermGlobalVarFlag )
				//   return Lv; // immediately, will not be used
				
				if (count > 10000) {
						Point pos = ray(ray.min_t);
						cout << "COUNT = " << count << " (Breaking) ray.min_t = " << ray.min_t << " max_t = "
						     << ray.max_t << " depth = " << ray.depth << " x = " << pos.x << " y = "
								 << pos.y << " z = " << pos.z << endl;
						if(count > 10005)
							break;
				}
#ifdef DEBUG
				p = ray(ray.min_t);
				cout << " VoronoiIntegrator::Li(iter=" << count << ") Lv.y = " << setw(5) << Lv.y()
						 << " Tr.y = " << setw(2) << Tr.y() << " ray.x = " << setw(5) << p.x 
						 << " ray.y = " << setw(5) << p.y << " ray.z = " << setw(5) << p.z << endl;
				cout << "  ( index = " << ray.index << " prev_index = " << ray.prev_index << " )" << endl;
#endif
				if (!scene->arepoMesh->AdvanceRayOneCellNew(ray, &t0, &t1, Lv, Tr, threadNum) )
						break;
				
        // roulette terminate ray marching if transmittance is small (only if not doing raw integrals)
        if (!Config.projColDens && Tr.y() < 1e-3) {
            const float continueProb = 0.5f;
						
#ifdef DEBUG
						Point p = ray(ray.min_t);
						cout << " roulette ray.x = " << p.x << " y = " << p.y << " z = " << p.z << " tr.y = " << Tr.y() << endl;
#endif
													
            if (rng.RandomFloat() > continueProb)
								break;
            Tr /= continueProb;
        }
		}
		
#ifdef DEBUG
		p = ray(ray.min_t);
		cout << " VoronoiIntegrator::Li(done_f) Lv.y = " << setw(6) << Lv.y()
				 << " Tr.y = " << Tr.y() << " ray.x = " << setw(5) << p.x 
				 << " ray.y = " << setw(5) << p.y << " ray.z = " << setw(5) << p.z << endl << endl;
#endif
    *T = Tr;
		return Lv;
}

VoronoiIntegrator *CreateVoronoiVolumeIntegrator()
{
    return new VoronoiIntegrator();
}

// ------------------------------- TreeSearchIntegrator -------------------------------
void TreeSearchIntegrator::Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer)
{
		// distribute rays to appropriate start tasks
		// based on tree topnodes
}

void TreeSearchIntegrator::RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene)
{
    tauSampleOffset = sample->Add1D(1);
    scatterSampleOffset = sample->Add1D(1);
}

Spectrum TreeSearchIntegrator::Transmittance(const Scene *scene, const Renderer *renderer, const Ray &ray,
																					   const Sample *sample, RNG &rng) const
{
		cout << "CHECK" << endl;
		endrun(1196);
		return Spectrum(1.0f); // added just to suppress return warning, CHECK
}

Spectrum TreeSearchIntegrator::Li(const Scene *scene, const Renderer *renderer, const Ray &ray,
															    const Sample *sample, RNG &rng, Spectrum *T,
																  int *prevEntryCell, int *prevEntryTetra, int threadNum) const
{
		IF_DEBUG(cout << "TreeSearchIntegrator::Li()" << endl);
		
    double t0, t1;
		
    if (!scene->arepoTree || !scene->arepoTree->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.0f) {
				IF_DEBUG(cout << " Returning! IntersectP t0 = " << t0 << " t1 = " << t1 << endl);
        *T = Spectrum(1.0f);
        return 0.0f;
    }
			
		// do emission only volume integration in AM
		Spectrum Lv(0.0);
		Spectrum Tr(1.0f);
		
		// propagate ray to arepo box, set exit point
		ray.min_t = t0;
		ray.max_t = t1;
		
#ifdef DEBUG
		ray(t0).print(" hitbox ");
		ray(t1).print(" exitbox ");
#endif		
		
		// maximum ray integration length from config
		if (Config.rayMaxT && Config.rayMaxT < ray.max_t)
				ray.max_t = Config.rayMaxT;
				
		IF_DEBUG(cout << " t0 = " << t0 << " t1 = " << t1
									<< " ray.min_t = " << ray.min_t << " ray.max_t = " << ray.max_t << endl);		
		

		// TODO: find initial tasks and exchange (or in preprocess)?
		ray.task = 0;
		ray.prevHSML = 10.0;
		
		// advance ray through box
#ifdef DEBUG
		int count = 0;

		Point p = ray(ray.min_t);
		cout << " TreeSearchIntegrator::Li(iter=" << count << ") Lv.y = " << setw(6) << Lv.y()
				 << " Tr.y = " << Tr.y() << " ray.x = " << setw(5) << p.x 
				 << " ray.y = " << setw(5) << p.y << " ray.z = " << setw(5) << p.z << endl;
#endif			 
		while( true ) {
				// TODO: SIGTERM
				// if( sigTermGlobalVarFlag )
				//   return Lv; // immediately, will not be used
#ifdef DEBUG
				count++;

				p = ray(ray.min_t);
				cout << " TreeSearchIntegrator::Li(iter=" << count << ") Lv.y = " << setw(5) << Lv.y()
						 << " Tr.y = " << setw(2) << Tr.y() << " ray.x = " << setw(5) << p.x 
						 << " ray.y = " << setw(5) << p.y << " ray.z = " << setw(5) << p.z << " (depth=" << ray.depth << ")" << endl;

				if (count > 10000) {
						cout << "COUNT > 10000 (Breaking)" << endl;
						break;
				}
#endif
				if (!scene->arepoTree->AdvanceRayOneStep(ray, &t0, &t1, Lv, Tr, threadNum) )
						break;
				
        // roulette terminate ray marching if transmittance is small (only if not doing raw integrals)
        if (!Config.projColDens && Tr.y() < 1e-3) {
            const float continueProb = 0.5f;
						
#ifdef DEBUG
						Point p = ray(ray.min_t);
						cout << " roulette ray.x = " << p.x << " y = " << p.y << " z = " << p.z << " tr.y = " << Tr.y() << endl;
#endif
													
            if (rng.RandomFloat() > continueProb)
								break;
            Tr /= continueProb;
        }
		}
#ifdef DEBUG
		p = ray(ray.min_t);
		cout << " TreeSearchIntegrator::Li(done_f) Lv.y = " << setw(6) << Lv.y()
				 << " Tr.y = " << Tr.y() << " ray.x = " << setw(5) << p.x 
				 << " ray.y = " << setw(5) << p.y << " ray.z = " << setw(5) << p.z << endl << endl;
#endif
    *T = Tr;
		return Lv;
}

TreeSearchIntegrator *CreateTreeSearchVolumeIntegrator()
{
    return new TreeSearchIntegrator();
}

// ------------------------------- QuadIntersectionIntegrator -------------------------------
void QuadIntersectionIntegrator::Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer)
{
		// calculate quads for hemisphere
		Point p0,p1,p2,p3;
		float theta[4], phi[4], rr[4];
		int N = 72; //36;
		
		for(int j=0; j < N/2; j++)
		{
			for(int i=-N/2; i < N/2; i++)
			{
				theta[0] = (j+0) * M_PI / N;
				theta[1] = (j+1) * M_PI / N;
				theta[2] = theta[1];
				theta[3] = theta[0];
				
				phi[0] = i * 2 * M_PI / N;
				phi[1] = phi[0];
				phi[2] = (i+1) * 2 * M_PI / N;
				phi[3] = phi[2];
										
				p0.x = sin(theta[0]) * sin(phi[0]);
				p0.y = sin(theta[0]) * cos(phi[0]);
				p0.z = cos(theta[0]);
				p1.x = sin(theta[1]) * sin(phi[1]);
				p1.y = sin(theta[1]) * cos(phi[1]);
				p1.z = cos(theta[1]);
				p2.x = sin(theta[2]) * sin(phi[2]);
				p2.y = sin(theta[2]) * cos(phi[2]);
				p2.z = cos(theta[2]);
				p3.x = sin(theta[3]) * sin(phi[3]);
				p3.y = sin(theta[3]) * cos(phi[3]);
				p3.z = cos(theta[3]);
				
				rr[0] = p0.x * p0.x + p0.y * p0.y;
				rr[1] = p1.x * p1.x + p1.y * p1.y;
				rr[2] = p2.x * p2.x + p2.y * p2.y;
				rr[3] = p3.x * p3.x + p3.y * p3.y;
				
				if (rr[0] > 1 || rr[1] > 1 || rr[2] > 1 || rr[3] > 1)
					continue;
					
				p0.z = sqrt(1.0 - rr[0]);
				p1.z = sqrt(1.0 - rr[1]);
				p2.z = sqrt(1.0 - rr[2]);
				p3.z = sqrt(1.0 - rr[3]);
				
				// save P0, compute S0 and S1 from 0-1 and 0-3, compute normal of S0xS1
				Vector s1_local = p0-p1;
				Vector s2_local = p1-p2;
				Vector qN_local = Normalize(Cross(s1_local,s2_local));
				
				P0.push_back(p2);
				S1.push_back(s1_local);
				S2.push_back(s2_local);
				QN.push_back(qN_local);
				
				nQuads++;
				
			} //i
		} //j
		
		// DEBUG:
		//nQuads = 1;
		//P0.push_back( Point(0.0,0.0,1.0) );
		//S1.push_back( Vector(0.0,0.5,0.0) );
		//S2.push_back( Vector(0.8,0.0,0.0) );
		//QN.push_back( Normalize(Cross(S1[0],S2[0])) );
				
		cout << "QuadIntersection::Preprocess(): made [" << nQuads << "] quads." << endl;
}

void QuadIntersectionIntegrator::RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene)
{
    tauSampleOffset = sample->Add1D(1);
    scatterSampleOffset = sample->Add1D(1);
}

Spectrum QuadIntersectionIntegrator::Transmittance(const Scene *scene, const Renderer *renderer, const Ray &ray,
																					   const Sample *sample, RNG &rng) const
{
		cout << "CHECK" << endl;
		endrun(1196);
		return Spectrum(1.0f); // added just to suppress return warning, CHECK
}

Spectrum QuadIntersectionIntegrator::Li(const Scene *scene, const Renderer *renderer, const Ray &ray,
															    const Sample *sample, RNG &rng, Spectrum *T,
																  int *prevEntryCell, int *prevEntryTetra, int threadNum) const
{
		IF_DEBUG(cout << "QuadIntersectionIntegrator::Li()" << endl);
		
		// config
		Spectrum Ledge = Spectrum::FromRGB(Config.rgbLine);
		Spectrum Lint  = Spectrum::FromRGB(Config.rgbTetra);
		Spectrum Lv(0.0);
		Spectrum Tr(1.0f);
		
		float edge_thickness = 0.005; // world coordinates
		float dn,a,d1_len,d2_len,d1_len_b,d2_len_b,I_edge;
		Vector d1,d2,P0P;
		Point Pi; // intersection point
		
		// loop over each quad
		for( int i=0; i < nQuads; i++ )
		{
			// check intersection with plane of quad
			dn = Dot(ray.d,QN[i]);
			
			if( dn > 0 )
				continue;
			
			// compute intersection coordinates of quad
			a = Dot(P0[i]-ray.o,QN[i]) / dn;
			Pi = ray.o + a * ray.d;
			
			P0P = Pi - P0[i];
			
			d1_len = Dot( P0P, S1[i] ) / S1[i].Length(); // scalar projection
			d2_len = Dot( P0P, S2[i] ) / S2[i].Length();
			
			// actually outside quad? then skip
			if( d1_len < 0 || d1_len > S1[i].Length() || d2_len < 0 || d2_len > S2[i].Length() )
			//if( d1_len < -edge_thickness || d1_len > S1[i].Length() || d2_len < 0 || d2_len > S2[i].Length() )
				continue;
			
			// no AA on outside edge, only inside edge, could get around by allowing the poly hit test 
			// to fail out to e.g. 3*sigma, then making sure the distances are the right sign below
			//d1_len = fabs(d1_len);
			//d2_len = fabs(d2_len);
			
			// compute distance to closest edge		
			d1_len_b = S1[i].Length() - d1_len;
			d2_len_b = S2[i].Length() - d2_len;
						
			d1_len = (d1_len < d2_len) ? d1_len : d2_len;
			d2_len = (d1_len_b < d2_len_b) ? d1_len_b : d2_len_b;
			d1_len = (d1_len < d2_len) ? d1_len : d2_len;
			
			// add to Lv as function of distance to nearest edge
			I_edge = exp(-2.0 * pow(d1_len/edge_thickness,2));
			
			Lv += Tr * I_edge * Ledge;
			Lv += Tr * (1.0-I_edge) * Lint; // interior of polygon

		}
			
    *T = Tr;
		return Lv;
}

QuadIntersectionIntegrator *CreateQuadIntersectionIntegrator()
{
    return new QuadIntersectionIntegrator();
}