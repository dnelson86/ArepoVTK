/*
 * renderer.cpp
 * dnelson
 */

#include <unistd.h>
 
#include "renderer.h"
#include "transform.h"
#include "arepo.h"
#include "camera.h"
#include "sampler.h"
#include "volume.h"
#include "integrator.h"

// RendererTask

void RendererTask::Run(int threadNum) {

    /* First Pass - Volume Ray Trace */
		Timer timer2;
		timer2.Start();
		
		// TODO: if taskFinishedArray[taskNum] == 1, this task is already done (from restart)
		//   immediate exit
		// TODO: if ( sigTermGlobalVarFlag == 1 ) then we are quickly cycling through the 
		//   remaining tasks, which we need to complete before we write the restart
		//   so: immediate exit
		
		//if( Config.filename.substr(0,4) == "none" && typeid(*renderer) != typeid(QuadIntersectionIntegrator) )
		//  return;
			
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
		float *rayWeights = new float[maxSamples];

    // Get samples from _Sampler_ and update image
    int sampleCount;
		
		// accelerate entry point location by using the entry point of the previous ray from this task
		int prevEntryCell  = -1;
		int prevEntryTetra = 0;
		
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0)
		{
				IF_DEBUG(cout << " [Thread=" << setw(2) << threadNum << " Task=" << setw(3) << taskNum 
				              << "] RendererTask::Run() maxSamples = " << maxSamples 
											<< " sampleCount = " << sampleCount << endl);
				
        // generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; i++)
				{
            // find camera ray for _sample[i]_
            rayWeights[i] = camera->GenerateRay(samples[i], &rays[i]);

						// find entry cell in mesh (done inside integrator for now - need t0,t1)
				}
				
				// Evaluate radiance along camera rays
				for (int i=0; i < sampleCount; i++)
				{
					if( rayWeights[i] > 0 )
						Ls[i] = rayWeights[i] * renderer->Li(scene, rays[i], &samples[i], rng, &Ts[i], 
						                                     &prevEntryCell, &prevEntryTetra, threadNum);
						// TODO:
						// if( sigTermGlobalVarFlag )
						//   break; (do not add any samples from this thread to the BlockArrays)
						
        }

        // Report sample results to _Sampler_, add contributions to image
        if (sampler->ReportResults(samples, rays, Ls, sampleCount))
        {
            for (int i = 0; i < sampleCount; ++i)
            {
							if( rayWeights[i] > 0 )
                camera->film->AddSample(samples[i], Ls[i], rays[i], threadNum);
            }
						// TODO:
						// mark taskFinishedArray[taskNum] = 1;
        }
    }

		float time = (float)timer2.Time();
		
		if( Config.verbose ) {
			cout << " [Thread=" << setw(2) << threadNum << " Task=" << setw(3) << taskNum 
					<< "] Raytracing phase: [" << setw(6) << time << "] seconds." << endl;
		}
		else {
			cout << ".";
			if( taskNum % 80 == 0 && taskCount > 100 )
				cout << "]" << endl << " [";
			fflush(stdout);
		}
		
    // Clean up after _SamplerRendererTask_ is done with its image region
    camera->film->UpdateDisplay(sampler->xPixelStart, sampler->yPixelStart, 
																sampler->xPixelEnd+1, sampler->yPixelEnd+1);
				
    delete sampler;
    delete[] samples;
    delete[] rays;
    delete[] Ls;
    delete[] Ts;
		delete[] rayWeights;
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

void Renderer::Render(const Scene *scene, int frameNum)
{
    volumeIntegrator->Preprocess(scene, camera, this);
    // Allocate and initialize _sample_
    Sample *sample = new Sample(sampler, volumeIntegrator, scene);
		
		// Make timer
		Timer *timer = new Timer();

    // image size and number of render tasks (threads) to run
    long long nPixels = camera->film->xResolution * camera->film->yResolution;
		int nTasks = Config.nTasks;
		long long nCores = numberOfCores();
		
		if (!nTasks) {
				// calculate number of tasks automatically
				int nTasks = max(TASK_MULT_FACT * nCores, nPixels / (TASK_MAX_PIXEL_SIZE*TASK_MAX_PIXEL_SIZE));
				Config.nTasks = RoundUpPowerOfTwo(nTasks);
		}
		
		cout << endl << "Render Setup: Using [" << Config.nTasks << "] tasks with [" << nCores 
				 << "] system cores (of " << sysconf(_SC_NPROCESSORS_ONLN) << " available) for " << camera->film->xResolution << "x"
		     << camera->film->yResolution << " image (" << nPixels << " pixels)." << endl << endl;
				 
		cout << "Rendering..." << endl << endl;
		timer->Start();
		
		// if not verbose, start single-line progress bar
		//if ( !Config.verbose )
		//	writeStatusBar(0, Config.nTasks);
		if( !Config.verbose )
			cout << " [";
		
		// add render tasks to queue
		vector<Task *> renderTasks;
		
		for (int i=0; i < Config.nTasks; i++)
				renderTasks.push_back(new RendererTask(scene, this, camera, sampler, sample, 
																							 nTasks-1-i, nTasks));

		// TODO: we are loading from a restart? if so load the taskFinishedArray now and fill it in
																							 
		// init tasks and run
		startTasks(renderTasks);
		
		// wait with this thread until workers are done
		waitUntilAllTasksDone();

		// free tasks when all done
		for (unsigned int i=0; i < renderTasks.size(); i++)
				delete renderTasks[i];
		
		if( !Config.verbose )
			cout << "]" << endl;
		
		// are we writing a restart and quitting due to SIGTERM?
		// TODO
		// if( sigTermGlobalVarFlag )
		//   write restart file: pixels/integrals BlockedArrays, taskFinishedArray, Config parameters (nTasks etc)
		//   return;
		
		// do rasterization stage with just one task
		if (Config.drawTetra || Config.drawVoronoi || Config.drawBBox || Config.drawSphere)
				Renderer::RasterizeStage(scene);		
		
    float seconds = (float)timer->Time();
		cout << endl << "[" << setw(3) << frameNum << "] Render complete (" 
		     << setw(6) << seconds << " seconds)." << endl << endl;
		
		// cleanup and store image, column integrals, and raw floats
		TasksCleanup();
    delete sample;
		
    camera->film->WriteImage(frameNum);
		camera->film->WriteIntegrals();
    IF_DEBUG(camera->film->WriteRawRGB());
}

void Renderer::RasterizeStage(const Scene *scene)
{
		Timer timer2;
		timer2.Start();
		
    /* Second Pass - Rasterize Mesh Faces */
    vector<Line> edges;
		
		// testing: rasterize edges of AM tetra
		if (Config.drawTetra && scene->arepoMesh) {
				Spectrum Ltetra	= Spectrum::FromRGB(Config.rgbTetra);
				for (int i = 0; i < scene->arepoMesh->Ndt; i++) {
						edges.clear();
						
						if(scene->arepoMesh->TetraEdges(i,&edges)) {
								for (unsigned int j = 0; j < edges.size(); j++) {
										IF_DEBUG(cout << " RL p1.x = " << edges[j].p1.x << " p1.y = " << edges[j].p1.y
																	<< " p1.z = " << edges[j].p1.z << " p2.x = " << edges[j].p2.x
																	<< " p2.y = " << edges[j].p2.y << " p2.z = " << edges[j].p2.z << endl);
										camera->RasterizeLine(edges[j].p1,edges[j].p2,Ltetra);
								}
						}
				} 
		} 

    // testing: rasterize Voronoi faces of AM
		if (Config.drawVoronoi && scene->arepoMesh) {
				edges.clear();
				Spectrum Lvor   = Spectrum::FromRGB(Config.rgbVoronoi);
				
				// compute edges
				int numFaces = scene->arepoMesh->ComputeVoronoiEdges();	
				
				for (int i = 0; i < numFaces-1; i++) {
						if (scene->arepoMesh->VoronoiEdges(i,&edges)) {
								for (unsigned int j = 0; j < edges.size(); j++) {
										IF_DEBUG(cout << " VE p1.x = " << edges[j].p1.x << " p1.y = " << edges[j].p1.y
																	<< " p1.z = " << edges[j].p1.z << " p2.x = " << edges[j].p2.x
																	<< " p2.y = " << edges[j].p2.y << " p2.z = " << edges[j].p2.z << endl);
										camera->RasterizeLine(edges[j].p1,edges[j].p2,Lvor);
								}

						}
				}
		} 

		// testing: acquire and rasterize edges of VR/AM BBox
		if (Config.drawBBox) {
				edges.clear();
				Spectrum Ledge  = Spectrum::FromRGB(Config.rgbLine);
				
				bool renderEdges = false;
				
				if( scene->arepoMesh )
					renderEdges = scene->arepoMesh->WorldBound().Edges(&edges);
				if( scene->arepoTree )
					renderEdges = scene->arepoTree->WorldBound().Edges(&edges);
					
				if (renderEdges) {
						for (unsigned int i = 0; i < edges.size(); ++i) {
								// rasterize individual line segment
								camera->RasterizeLine(edges[i].p1,edges[i].p2,Ledge);
						}
				}
		} 
		
		// draw test sphere
		if (Config.drawSphere) {
				edges.clear();
				Spectrum Ledge = Spectrum::FromRGB(Config.rgbLine);
				
				float rgb_r[3]; rgb_r[0] = 1.0; rgb_r[1] = 0.0; rgb_r[2] = 0.0;
				float rgb_g[3]; rgb_g[0] = 0.0; rgb_g[1] = 1.0; rgb_g[2] = 0.0;
				float rgb_b[3]; rgb_b[0] = 0.0; rgb_b[1] = 0.0; rgb_b[2] = 1.0;
				
				Spectrum Lr = Spectrum::FromRGB(rgb_r);
				Spectrum Lg = Spectrum::FromRGB(rgb_g);
				Spectrum Lb = Spectrum::FromRGB(rgb_b);
				
				Point p0,p1,p2,p3;
				float theta[4], phi[4], rr[4];
				
				// inspired by http://paulbourke.net/dome/fisheye/source1.c (Paul Bourke)
				int N = 36;
				
				for(int j=0; j < N/2; j++)
				{
					for(int i=-N/2; i < N/2; i++)
					{
						// OPTION A: (also switch j loop, start at -N_lat/2
						/*p0.y = j * 2 * N;
						p1.y = (j+1) * 2 * N;
						p2.y = p1.y;
						p3.y = p0.y;
						
						p0.x = i * 2 * N;
						p1.x = p0.x;
						p2.x = (i+1) * 2 * N;
						p3.x = p2.x;*/
						
						// OPTION B:
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
						
						// BOTH:						
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
						
						// draw the polygon with vertices p0,p1,p2,p3
						camera->RasterizeLine(p0,p1,Lr);
						camera->RasterizeLine(p1,p2,Lg);
						camera->RasterizeLine(p2,p3,Lb);
						camera->RasterizeLine(p3,p0,Ledge);
						
					} //i
				} //j
		}

		float time = (float)timer2.Time();
		cout << endl << " [Task=00] Rasterization phase: [" << time << "] seconds." << endl;
}

Spectrum Renderer::Li(const Scene *scene, const Ray &ray, const Sample *sample, RNG &rng, Spectrum *T, int *prevEntryCell, int *prevEntryTetra, int taskNum) const
{
    Spectrum Lvi = volumeIntegrator->Li(scene, this, ray, sample, rng, T, prevEntryCell, prevEntryTetra, taskNum);
		
    return Lvi;
}

Spectrum Renderer::Transmittance(const Scene *scene, const Ray &ray, const Sample *sample, RNG &rng) const
{
    return volumeIntegrator->Transmittance(scene, this, ray, sample, rng);
}
