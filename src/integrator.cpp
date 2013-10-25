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

// EmissionIntegrator
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

// VoronoiIntegrator
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
		ray.min_t = t0 /*+ 0.0*/;
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
#ifdef DEBUG
		int count = 0;

		Point p = ray(ray.min_t);
		cout << " VoronoiIntegrator::Li(iter=" << count << ") Lv.y = " << setw(6) << Lv.y()
				 << " Tr.y = " << Tr.y() << " ray.x = " << setw(5) << p.x 
				 << " ray.y = " << setw(5) << p.y << " ray.z = " << setw(5) << p.z << endl;
#endif			 
		while( true ) {
#ifdef DEBUG
				count++;

				p = ray(ray.min_t);
				cout << " VoronoiIntegrator::Li(iter=" << count << ") Lv.y = " << setw(5) << Lv.y()
						 << " Tr.y = " << setw(2) << Tr.y() << " ray.x = " << setw(5) << p.x 
						 << " ray.y = " << setw(5) << p.y << " ray.z = " << setw(5) << p.z << endl;
				cout << "  ( index = " << ray.index << " prev_index = " << ray.prev_index << " )" << endl;

				if (count > 10000) {
						cout << "COUNT > 10000 (Breaking)" << endl;
						break;
				}
#endif
				if (!scene->arepoMesh->AdvanceRayOneCellNew(ray, &t0, &t1, Lv, Tr, threadNum) )
						break;
				
        // roulette terminate ray marching if transmittance is small
        if (Tr.y() < 1e-3) {
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
