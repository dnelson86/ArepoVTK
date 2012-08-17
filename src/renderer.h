/*
 * renderer.h
 * dnelson
 */
 
#ifndef AREPO_RT_RENDERER_H
#define AREPO_RT_RENDERER_H

#include "ArepoRT.h"
#include "util.h"

// Renderer
class Renderer {
public:
    // construction
    Renderer(Sampler *s, Camera *c, VolumeIntegrator *vi);
    ~Renderer();
		
		// methods
    void Render(const Scene *scene);
		void RasterizeStage(const Scene *scene);
		
    Spectrum Li(const Scene *scene, const Ray &ray, const Sample *sample, RNG &rng, Spectrum *T = NULL, int prevEntryCell = -1) const;
    Spectrum Transmittance(const Scene *scene, const Ray &ray, const Sample *sample, RNG &rng) const;
		
private:
    // data
    Sampler *sampler;
    Camera *camera;
    VolumeIntegrator *volumeIntegrator;
};

// RendererTask
class RendererTask : public Task {
public:
    // construction
    RendererTask(const Scene *sc, Renderer *ren, Camera *c, Sampler *ms, Sample *sam, int tn, int tc)
    {
				IF_DEBUG(cout << "RendererTask(" << tn << ", " << tc << ") constructor." << endl);

        scene = sc; renderer = ren; camera = c; mainSampler = ms;
        origSample = sam; taskNum = tn; taskCount = tc;
				prevEntryCell = -1;
    }
    void Run();
		
private:
    // data
    const Renderer *renderer;
		const Scene *scene;
    Camera *camera;
    Sampler *mainSampler;
    Sample *origSample;
    int taskNum, taskCount;
		int prevEntryCell;
};



#endif
