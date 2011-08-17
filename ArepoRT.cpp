/*
 * ArepoRT.cpp
 * dnelson
 */
 
#include "ArepoRT.h"
#include "integrator.h"
#include "renderer.h"
#include "volume.h"
#include "arepo.h"

void rtTestRenderScene(string filename)
{
		// config - camera
		Point cameraPos    = Point(10.5,10.5,0.0); //on z=0 plane
		Point cameraLook   = Point(10.5,10.5,1.0);
		//Point cameraPos    = Point(14.0,11.6,2.0);
		//Point cameraLook   = Point(10.5,10.5,3.5); //center of box
		Vector cameraVecUp = Vector(0.0,1.0,0.0);
		
		// set transform
		Transform volume2world;
		Transform world2camera;

		volume2world = volume2world * Translate(Vector(10.0,10.0,40.0));
		world2camera = world2camera * LookAt(cameraPos, cameraLook, cameraVecUp);
		
		IF_DEBUG(world2camera.print("world2camera:"));
		IF_DEBUG(Inverse(world2camera).print("camera2world:"));
		IF_DEBUG(volume2world.print("volume2world:"));
		
		// MakeCamera(): filter and film
		Filter *filter       = CreateBoxFilter();
		Film *film           = CreateFilm(filter);
		
		// ::MakeCamera(): camera
		Camera *camera       = CreateOrthographicCamera(Inverse(world2camera), film);
		
		// Sampler
		Sampler *sampler     = CreateStratifiedSampler(film, camera);
		
		// Integrator
		//VolumeIntegrator *vi = CreateEmissionVolumeIntegrator(Config.viStepSize);
		VolumeIntegrator *vi  = CreateVoronoiVolumeIntegrator();
		
		// renderer
		Renderer *re         = new Renderer(sampler, camera, vi);
		
		// create volume/density/scene geometry
		//VolumeRegion *vr     = CreateGridVolumeRegion(volume2world, filename);
		VolumeRegion *vr      = NULL;
		
		// voronoi mesh
		ArepoMesh *arepoMesh  = CreateArepoMesh(volume2world);
		//arepoMesh->DumpMesh();

		// scene
		Scene *scene          = new Scene(vr,arepoMesh);
		
#ifdef DEBUG
    if(arepoMesh) arepoMesh->WorldBound().print("ArepoMesh WorldBound ");
		if(vr)        vr->WorldBound().print("VolumeRegion WorldBound ");
		
		//Point testp(0.2, 0.2, 0.2);
		//cout << "VR Test @ (" << testp.x << " " << testp.y << " " << testp.z << ") result = "
		//		 << vr->Density(testp) << endl;
#endif
				 
		// render
		if (filter && film && camera && sampler && vi && scene && re)
				re->Render(scene);
				

		delete re;
		delete scene;
}

/*
    Renderer *renderer = renderOptions->MakeRenderer();
		  Renderer *renderer = NULL;
			Camera *camera = MakeCamera();
			  Filter *filter = MakeFilter(FilterName, FilterParams);
				Film *film = MakeFilm(FilmName, FilmParams, filter);
				Camera *camera = ::MakeCamera(CameraName, CameraParams,CameraToWorld, renderOptions->transformStartTime,
				                              renderOptions->transformEndTime, film);
			Sampler *sampler = MakeSampler(SamplerName, SamplerParams, camera->film, camera);
			VolumeIntegrator *volumeIntegrator = MakeVolumeIntegrator(VolIntegratorName,VolIntegratorParams);
			renderer = new SamplerRenderer(sampler, camera, surfaceIntegrator,volumeIntegrator, visIds);
    Scene *scene = renderOptions->MakeScene();
		  VolumeRegion *volumeRegion;
			Scene *scene = new Scene(accelerator, lights, volumeRegion);
    if (scene && renderer) renderer->Render(scene);
    TasksCleanup();
    delete renderer;
    delete scene;
*/

ConfigStruct Config;

int main (int argc, char* argv[])
{
		string filename;

		cout << "ArepoRT v" << AREPO_RT_VERSION << " (" << __DATE__ << ")\n" 
				 << "--------------------------\n";
		IF_DEBUG(cout << "DEBUGGING ENALBED.\n\n");
		
		// command line arguments
	/*	if (argc != 2) {
				cout << "usage: ArepoRT <scenefile.txt>\n";
				return 0;
		} else {
				Config.filename = argv[1];
		} */

		// read config
		
		// init
#ifdef ENABLE_AREPO
    cout << "AREPO ENABLED." << endl;

    Arepo arepo = Arepo(Config.filename, Config.paramFilename);
		
    arepo.Init(&argc,&argv);

    arepo.LoadSnapshot();
#endif		

		// debug test render
		rtTestRenderScene(Config.filename);
		
		// cleanup
#ifdef ENABLE_AREPO
		arepo.Cleanup();
#endif

		return 0;
}
