/*
 * integrator.h
 * dnelson
 */
 
#ifndef AREPO_RT_INTEGRATOR_H
#define AREPO_RT_INTEGRATOR_H

#include "ArepoRT.h"

class VolumeIntegrator {
public:
		// pure virtual
    virtual ~VolumeIntegrator() { };
		
    virtual void Preprocess(const Scene *scene, const Camera *camera,const Renderer *renderer) {
    }
    virtual void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene) {
    }
		
    virtual Spectrum Li(const Scene *scene, const Renderer *renderer, const Ray &ray, 
		                    const Sample *sample, RNG &rng, Spectrum *transmittance, 
												int *prevEntryCell, int *prevEntryTetra, int taskNum) const = 0;
    virtual Spectrum Transmittance(const Scene *scene, const Renderer *renderer, const Ray &ray,
																	 const Sample *sample, RNG &rng) const = 0;
};

class EmissionIntegrator : public VolumeIntegrator {
public:
    // construction
    EmissionIntegrator(float ss) {
				IF_DEBUG(cout << "EmissionIntegrator(" << ss << ") constructor." << endl);
				stepSize = ss;
    }
		//~EmissionIntegrator() { };

    // methods
		void Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer) { }
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    Spectrum Li(const Scene *scene, const Renderer *renderer, const Ray &ray, 
		            const Sample *sample, RNG &rng, Spectrum *transmittance, int *prevEntryCell, int *prevEntryTetra, int taskNum) const;
    Spectrum Transmittance(const Scene *scene, const Renderer *,
            const Ray &ray, const Sample *sample, RNG &rng) const;
private:
    // data
    float stepSize;
    int tauSampleOffset, scatterSampleOffset;
};

class VoronoiIntegrator : public VolumeIntegrator {
public:
    // consturction
    VoronoiIntegrator() {
				IF_DEBUG(cout << "VoronoiIntegrator() constructor." << endl);
    }
		//~VoronoiIntegrator() { };

    // methods
		void Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer);
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    Spectrum Li(const Scene *scene, const Renderer *renderer, const Ray &ray, 
		            const Sample *sample, RNG &rng, Spectrum *transmittance, int *prevEntryCell, int *prevEntryTetra, int taskNum) const;
    Spectrum Transmittance(const Scene *scene, const Renderer *,
            const Ray &ray, const Sample *sample, RNG &rng) const;
private:
    // data
    int tauSampleOffset, scatterSampleOffset;
};

class TreeSearchIntegrator : public VolumeIntegrator {
public:
    // consturction
    TreeSearchIntegrator() {
				IF_DEBUG(cout << "TreeSearchIntegrator() constructor." << endl);
    }
		//~TreeSearchIntegrator() { };

    // methods
		void Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer);
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    Spectrum Li(const Scene *scene, const Renderer *renderer, const Ray &ray, 
		            const Sample *sample, RNG &rng, Spectrum *transmittance, int *prevEntryCell, int *prevEntryTetra, int threadNum) const;
    Spectrum Transmittance(const Scene *scene, const Renderer *,
            const Ray &ray, const Sample *sample, RNG &rng) const;
private:
    // data
    int tauSampleOffset, scatterSampleOffset;
};

class QuadIntersectionIntegrator : public VolumeIntegrator {
public:
    // consturction
    QuadIntersectionIntegrator() {
				IF_DEBUG(cout << "QuadIntersectionIntegrator() constructor." << endl);
    }
		//~QuadIntersectionIntegrator() { };

    // methods
		void Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer);
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    Spectrum Li(const Scene *scene, const Renderer *renderer, const Ray &ray, 
		            const Sample *sample, RNG &rng, Spectrum *transmittance, int *prevEntryCell, int *prevEntryTetra, int threadNum) const;
    Spectrum Transmittance(const Scene *scene, const Renderer *,
            const Ray &ray, const Sample *sample, RNG &rng) const;
private:
    // data
    int tauSampleOffset, scatterSampleOffset;
		int nQuads;
		vector<Point> P0;
		vector<Vector> S1;
		vector<Vector> S2;
		vector<Vector> QN;
};

EmissionIntegrator *CreateEmissionVolumeIntegrator(const float stepSize);
VoronoiIntegrator *CreateVoronoiVolumeIntegrator();
TreeSearchIntegrator *CreateTreeSearchVolumeIntegrator();
QuadIntersectionIntegrator *CreateQuadIntersectionIntegrator();

#endif
