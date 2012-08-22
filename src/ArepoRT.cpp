/*
 * ArepoRT.cpp
 * dnelson
 */
 
#include "ArepoRT.h"

#include "sampler.h"
#include "transform.h"
#include "camera.h"

#include "voronoi_3db.h"
#include "integrator.h"
#include "volume.h"
#include "renderer.h"

#include "transfer.h"
#include "util.h"
#include "spectrum.h"

void rtIsoDiskCosmoCutoutRender()
{
	Timer timer;
	timer.Start();
	
	const Spectrum sig_a = Spectrum::FromRGB(Config.rgbAbsorb);

	// constant setup		
	VolumeIntegrator *vi  = CreateVoronoiVolumeIntegrator();	
	TransferFunction *tf  = new TransferFunction(sig_a);
	ArepoMesh *arepoMesh  = new ArepoMesh(tf);
	Scene *scene          = new Scene(NULL, arepoMesh);
	
	// setup camera
	//Point cameraPos    = Point(boxSize/2.0,boxSize/2.0,0.0); //face on		
	//Point cameraPos    = Point(boxSize,boxSize/2.0,boxSize/2.0); //edge on (#1)
	//Point cameraPos    = Point(0.0,boxSize/2.0,boxSize/2.0); //edge on (#2)
	//Point cameraLook   = Point(boxSize/2.0,boxSize/2.0,boxSize/2.0); //center of box
	//Point cameraLook   = cameraPos+Point(-1.0,0.0,0.0); //20mpc debug
	
	Vector cameraVecUp     = Vector(0.0,1.0,0.0);
	Transform world2camera = LookAt(Point(Config.cameraPosition), Point(Config.cameraLookAt), cameraVecUp);
	
	// setup transfer function
	for (unsigned int i=0; i < Config.tfSet.size(); i++)
		tf->AddParseString(Config.tfSet[i]);

	// would recreate per frame
	Filter *filter       = CreateBoxFilter();
	Film *film           = CreateFilm(filter);
	Camera *camera       = NULL;
	
	if (Config.cameraFOV)	camera = CreatePerspectiveCamera(Inverse(world2camera), film);
	else								  camera = CreateOrthographicCamera(Inverse(world2camera), film);
		
	Sampler *sampler     = CreateStratifiedSampler(film, camera);
	Renderer *re         = new Renderer(sampler, camera, vi);
	
	cout << endl << "Raytracer Init and Arepo Preprocessing: [" << (float)timer.Time() << "] seconds." << endl;
	
	// render
	if (re && scene)
		re->Render(scene);
			
	delete re;
	delete scene;
}

void rtTestRenderScene(string filename)
{
	// config - camera
	Vector cameraVecUp     = Vector(0.0,1.0,0.0);
	Transform world2camera = LookAt(Point(Config.cameraPosition), Point(Config.cameraLookAt), cameraVecUp);
	
	IF_DEBUG(world2camera.print("world2camera:"));
	IF_DEBUG(Inverse(world2camera).print("camera2world:"));
	
	// Camera (Filter, Film), Sampler
	Filter *filter       = CreateBoxFilter();
	Film *film           = CreateFilm(filter);
	Camera *camera       = NULL;
	
	if (Config.cameraFOV)
		camera = CreatePerspectiveCamera(Inverse(world2camera), film);
	else
		camera = CreateOrthographicCamera(Inverse(world2camera), film);
		
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
	Spectrum s1 = Spectrum::FromNamed("red");
	Spectrum s2 = Spectrum::FromNamed("green");
	Spectrum s3 = Spectrum::FromNamed("blue");
	tf->AddConstant(TF_VAL_DENS,s3);
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
	for (int i=0; i < NumGas; i++) {
		SphP[i].Density      = 0.01;
		SphP[i].Grad.drho[0] = 0.0;
		SphP[i].Grad.drho[1] = 0.0;
		SphP[i].Grad.drho[2] = 0.0;
	}
	
	// Arepo2b overrides
	SphP[6].Density      = 20.0;  // center point //3.0
	SphP[6].Grad.drho[0] = 10.0;
	SphP[6].Grad.drho[1] = 10.0;
	SphP[6].Grad.drho[2] = 2.0;
	
	//SphP[2].Density      = 20.0; // upper right near corner
	//SphP[5].Density      = 2.8; // upper right far corner
	
#ifdef DEBUG
	for (int i=0; i < NumGas; i++) {
			cout << "SphP[" << setw(2) << i << "] dens = " << setw(10) << SphP[i].Density 
			 << " grad.x = " << setw(10) << SphP[i].Grad.drho[0] << " grad.y = " 
			 << setw(10) << SphP[i].Grad.drho[1] << " grad.z = " << setw(10) << SphP[i].Grad.drho[2] << endl;
	}
	
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
	//rtTestRenderScene(Config.filename);
	rtIsoDiskCosmoCutoutRender();
	
	// cleanup
#ifdef ENABLE_AREPO
	arepo.Cleanup();
#endif

	return 0;
}
