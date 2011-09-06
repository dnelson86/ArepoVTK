/*
 * ArepoRT.cpp
 * dnelson
 */
 
#include "ArepoRT.h"
#include "integrator.h"
#include "renderer.h"
#include "volume.h"
#include "arepo.h"
#include "transfer.h"

void rtTestIsoDiskRender()
{
		const Spectrum sig_a = Spectrum::FromRGB(Config.rgbAbsorb);

		// constant setup		
		VolumeIntegrator *vi  = CreateVoronoiVolumeIntegrator();	
		TransferFunction *tf  = new TransferFunction(sig_a);
		ArepoMesh *arepoMesh  = new ArepoMesh(tf);
		Scene *scene          = new Scene(NULL, arepoMesh);
		
		// setup camera
		Point cameraPos    = Point(boxSize/2.0,boxSize/2.0,5000.0); //face on
		//Point cameraPos      = Point(11700.0,11200.0,0.0); //20mpc debug 1
		//Point cameraPos      = Point(400.0,15000.0,0.0); //20mpc debug 2
		
		//Point cameraPos    = Point(boxSize,boxSize/2.0,boxSize/2.0); //edge on
		//Point cameraPos    = Point(boxSize/2.0,0.0,boxSize*0.25); //high i angled
		
		//Point cameraLook   = Point(boxSize/2.0,boxSize/2.0,boxSize/2.0); //center of box
		Point cameraLook = cameraPos+Point(0.0,0.0,1.0); //20mpc debug
		Vector cameraVecUp = Vector(0.0,1.0,0.0);

		Transform world2camera = LookAt(cameraPos, cameraLook, cameraVecUp);
		
		// setup transfer function
		for (int i=0; i < Config.tfSet.size(); i++)
				tf->AddParseString(Config.tfSet[i]);

		// would recreate per frame
		Filter *filter       = CreateBoxFilter();
		Film *film           = CreateFilm(filter);
		Camera *camera       = CreateOrthographicCamera(Inverse(world2camera), film);
		Sampler *sampler     = CreateStratifiedSampler(film, camera);
		Renderer *re         = new Renderer(sampler, camera, vi);
		
		// render
		if (re && scene)
				re->Render(scene);
				
		delete re;
		delete scene;
}

void rtTestRenderScene(string filename)
{
		// config - camera
		//Point cameraPos    = Point(0.5,0.5,0.0); //centered on z=0 plane edge, x/y axis aligned
		//Point cameraPos      = Point(1.0,0.5,0.5); //centered on x=1.0 plane edge
		//Point cameraPos    = Point(2.0,2.0,-2.0); //angled above
		
		Point cameraPos = Point(0.4,0.4,0.0);
		
		//Point cameraLook   = Point(0.5,0.5,0.5); //center of box
		Point cameraLook = Point(0.4,0.4,1.0);
		Vector cameraVecUp = Vector(0.0,1.0,0.0);
		
		// set transform
		Transform world2camera = LookAt(cameraPos, cameraLook, cameraVecUp);
		
		IF_DEBUG(world2camera.print("world2camera:"));
		IF_DEBUG(Inverse(world2camera).print("camera2world:"));
		
		// Camera (Filter, Film), Sampler
		Filter *filter       = CreateBoxFilter();
		Film *film           = CreateFilm(filter);
		Camera *camera       = CreateOrthographicCamera(Inverse(world2camera), film);
		Sampler *sampler     = CreateStratifiedSampler(film, camera);
		
		// Integrator
		//VolumeIntegrator *vi = CreateEmissionVolumeIntegrator(Config.viStepSize);
		VolumeIntegrator *vi  = CreateVoronoiVolumeIntegrator();
		
		// renderer
		Renderer *re         = new Renderer(sampler, camera, vi);
		
		// transfer function
		const Spectrum sig_a = Spectrum::FromRGB(Config.rgbAbsorb);
		TransferFunction *tf = new TransferFunction(sig_a);
		
		// debugging:
		Spectrum s1 = Spectrum::FromRGB(Config.rgbEmit);
		Spectrum s2 = Spectrum::FromNamed("green");
		Spectrum s3 = Spectrum::FromNamed("blue");
		//tf->AddConstant(TF_VAL_DENS,s3);
		//tf->AddTophat(TF_VAL_DENS,5.0,10.0,spec);
		tf->AddGaussian(TF_VAL_DENS,2.8,0.1,s1);
		tf->AddGaussian(TF_VAL_DENS,5.0,0.5,s2);
		
		// create volume/density/scene geometry (debugging only)
		//VolumeRegion *vr     = CreateGridVolumeRegion(volume2world, filename);
		VolumeRegion *vr      = NULL;
		
		// voronoi mesh
		ArepoMesh *arepoMesh  = new ArepoMesh(tf);

		// scene
		Scene *scene          = new Scene(vr, arepoMesh);
				
		// debugging only (Arepo2b overrides)
		for (int i=0; i < N_gas; i++) {
				SphP[i].Density      = 0.01;
				SphP[i].Grad.drho[0] = 0.0;
				SphP[i].Grad.drho[1] = 0.0;
				SphP[i].Grad.drho[2] = 0.0;
		}
		
		// Arepo2b overrides
		SphP[6].Density      = 3.0;  // center point
		SphP[6].Grad.drho[0] = 10.0;
		SphP[6].Grad.drho[1] = 10.0;
		SphP[6].Grad.drho[2] = 2.0;
		
		SphP[2].Density      = 20.0; // upper right near corner
		
		SphP[5].Density      = 2.8; // upper right far corner
		
		for (int i=0; i < N_gas; i++) {
				cout << "SphP[" << setw(2) << i << "] dens = " << setw(10) << SphP[i].Density 
						 << " grad.x = " << setw(10) << SphP[i].Grad.drho[0] << " grad.y = " 
						 << setw(10) << SphP[i].Grad.drho[1] << " grad.z = " << setw(10) << SphP[i].Grad.drho[2] << endl;
	  }
		
#ifdef DEBUG
    if(arepoMesh) arepoMesh->WorldBound().print("ArepoMesh WorldBound ");
		if(vr)        vr->WorldBound().print("VolumeRegion WorldBound ");
#endif
				 
		// render
		if (scene && re)
				re->Render(scene);
				
		delete re;
		delete scene;
}

ConfigSet Config;

int main (int argc, char* argv[])
{
		cout << "ArepoRT v" << AREPO_RT_VERSION << " (" << __DATE__ << ")\n" 
				 << "--------------------------\n";
		IF_DEBUG(cout << "DEBUGGING ENALBED.\n\n");
		
		// command line arguments
		if (argc != 2) {
				cout << "usage: ArepoRT <configfile.txt>\n";
				return 0;
		} else {
				// read config
				Config.ReadFile( argv[1] );
				if (Config.verbose)
						Config.print();
		} 
		
		// init
#ifdef ENABLE_AREPO
    Arepo arepo = Arepo(Config.filename, Config.paramFilename);
		
    arepo.Init(&argc,&argv);
    arepo.LoadSnapshot();
#endif

		// debug test render
		rtTestRenderScene(Config.filename);
		
		//rtTestIsoDiskRender();
		
		// cleanup
#ifdef ENABLE_AREPO
		arepo.Cleanup();
#endif

		return 0;
}
