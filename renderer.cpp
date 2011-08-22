/*
 * renderer.cpp
 * dnelson
 */
 
#include "renderer.h"

// RendererTask

void RendererTask::Run() {

    /* First Pass - Volume Ray Trace */
		Timer timer2;
		timer2.Start();
		
    // Get sub-_Sampler_ for _SamplerRendererTask_
    Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
    if (!sampler) {
        return;
    }

    // Declare local variables used for rendering loop
    RNG rng(taskNum);

    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);
    Ray *rays = new Ray[maxSamples];
    Spectrum *Ls = new Spectrum[maxSamples];
    Spectrum *Ts = new Spectrum[maxSamples];

    // Get samples from _Sampler_ and update image
    int sampleCount;
		
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0)
		{
				IF_DEBUG(cout << " RendererTask::Run() maxSamples = " << maxSamples 
											<< " sampleCount = " << sampleCount << endl);
				
        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i)
				{
            // Find camera ray for _sample[i]_
            float rayWeight = camera->GenerateRay(samples[i], &rays[i]);

						// Evaluate radiance along camera ray
						if (rayWeight > 0.0f)
								Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng, &Ts[i]);
						else {
								Ls[i] = 0.0f;
								Ts[i] = 1.0f;
						}

						// error check on value of Ls[i]
        }

        // Report sample results to _Sampler_, add contributions to image
        if (sampler->ReportResults(samples, rays, Ls, sampleCount))
        {
            for (int i = 0; i < sampleCount; ++i)
            {
                camera->film->AddSample(samples[i], Ls[i]);
            }
        }
    }

		float time = (float)timer2.Time();
		cout << " Raytracing phase: [" << time << "] seconds." << endl;
		timer2.Reset();
		timer2.Start();
		
    /* Second Pass - Rasterize Mesh Faces */
    vector<Line> edges;
		
		// testing: rasterize edges of AM tetra
		if (Config.drawTetra) {
				Spectrum Ltetra	= Spectrum::FromRGB(Config.rgbTetra);
				for (int i = 0; i < scene->arepoMesh->Ndt; i++) {
						edges.clear();
						
						if(scene->arepoMesh->TetraEdges(i,&edges)) {
								for (int j = 0; j < edges.size(); j++) {
										IF_DEBUG(cout << " RL p1.x = " << edges[j].p1.x << " p1.y = " << edges[j].p1.y
																	<< " p1.z = " << edges[j].p1.z << " p2.x = " << edges[j].p2.x
																	<< " p2.y = " << edges[j].p2.y << " p2.z = " << edges[j].p2.z << endl);
										camera->RasterizeLine(edges[j].p1,edges[j].p2,Ltetra);
								}
						}
				} 
		} 

    // testing: rasterize Voronoi faces of AM
		if (Config.drawVoronoi) {
				edges.clear();
				Spectrum Lvor   = Spectrum::FromRGB(Config.rgbVoronoi);
				
				//for (int i = 0; i < scene->arepoMesh->Nvf; i++) {
				int i = 26;
						if (scene->arepoMesh->VoronoiEdges(i,&edges)) {
								for (int j = 0; j < edges.size(); j++) {
										cout << " VE p1.x = " << edges[j].p1.x << " p1.y = " << edges[j].p1.y
																	<< " p1.z = " << edges[j].p1.z << " p2.x = " << edges[j].p2.x
																	<< " p2.y = " << edges[j].p2.y << " p2.z = " << edges[j].p2.z << endl;
										camera->RasterizeLine(edges[j].p1,edges[j].p2,Lvor);
								}

						}
				//}
		} 

		// testing: acquire and rasterize edges of VR/AM BBox
		if (Config.drawBBox) {
				edges.clear();
				Spectrum Ledge  = Spectrum::FromRGB(Config.rgbLine);
				
				if (scene->arepoMesh->WorldBound().Edges(&edges)) {
						for (int i = 0; i < edges.size(); ++i) {
								// rasterize individual line segment
								camera->RasterizeLine(edges[i].p1,edges[i].p2,Ledge);
						}
				}
		} 

		time = (float)timer2.Time();
		cout << " Rasterization phase: [" << time << "] seconds." << endl;		
		
    // Clean up after _SamplerRendererTask_ is done with its image region
    camera->film->UpdateDisplay(sampler->xPixelStart, sampler->yPixelStart, 
																sampler->xPixelEnd+1, sampler->yPixelEnd+1);
				
    delete sampler;
    delete[] samples;
    delete[] rays;
    delete[] Ls;
    delete[] Ts;
}

// Renderer

Renderer::Renderer(Sampler *s, Camera *c, VolumeIntegrator *vi)
{
		IF_DEBUG(cout << "Renderer() constructor." << endl);

    sampler = s;
    camera = c;
    volumeIntegrator = vi;
}

Renderer::~Renderer()
{
    delete sampler;
    delete camera;
    delete volumeIntegrator;
}

void Renderer::Render(const Scene *scene)
{
    volumeIntegrator->Preprocess(scene, camera, this);
    // Allocate and initialize _sample_
    Sample *sample = new Sample(sampler, volumeIntegrator, scene);
		
		// Make timer
		Timer *timer = new Timer();

    // Create and launch _SamplerRendererTask_s for rendering image
    // Compute number of _SamplerRendererTask_s to create for rendering
    int nPixels = camera->film->xResolution * camera->film->yResolution;
    //int nTasks = max(32 * NumSystemCores(), nPixels / (16*16));
    //nTasks = RoundUpPow2(nTasks);
		int nTasks = Config.nTasks;
		int i=0;
		
		cout << endl << "Render Setup: Using [" << nTasks << "] tasks for " << camera->film->xResolution << "x"
		     << camera->film->yResolution << " image (" << nPixels << " pixels)." << endl << endl;
				 
		cout << "Rendering..." << endl << endl;
		
		timer->Start();
		
    RendererTask *renderTask = new RendererTask(scene, this, camera, sampler, sample, nTasks-1-i, nTasks);

    renderTask->Run();
		
    float seconds = (float)timer->Time();
		
		cout << endl << "Render complete (" << seconds << " seconds)." << endl << endl;
		
		// cleanup and store image
    delete renderTask;
    delete sample;
		
    camera->film->WriteImage();
		IF_DEBUG(camera->film->WriteRawRGB());
}

Spectrum Renderer::Li(const Scene *scene, const Ray &ray, const Sample *sample, RNG &rng, Spectrum *T) const
{
    Spectrum Lvi = volumeIntegrator->Li(scene, this, ray, sample, rng, T);
		
    return Lvi;
}

Spectrum Renderer::Transmittance(const Scene *scene, const Ray &ray, const Sample *sample, RNG &rng) const
{
    return volumeIntegrator->Transmittance(scene, this, ray, sample, rng);
}
