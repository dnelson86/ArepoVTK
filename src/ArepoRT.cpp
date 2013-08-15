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
#include "keyframe.h"

void rtIsoDiskCosmoCutoutRender()
{
	Timer timer;
	timer.Start();
	
	const Spectrum sig_a = Spectrum::FromRGB(Config.rgbAbsorb);

	// constant setup		
	TransferFunction *tf  = new TransferFunction(sig_a);
	ArepoMesh *arepoMesh  = new ArepoMesh(tf);
	Scene *scene          = new Scene(NULL, arepoMesh);
	FrameManager *fm      = new FrameManager(Config.kfSet);
	
	// setup transfer function
	for (unsigned int i=0; i < Config.tfSet.size(); i++)
		tf->AddParseString(Config.tfSet[i]);

#ifdef SPECIAL_BOUNDARY
	// override density of ID=-2 (inner spoon boundary cells) with density much higher than surrounding gas
	// put a TF near this value but below to highlight spoon boundary
	for(int i = 0; i < NumGas; i++)
	{
		if( P[i].ID == -2 )
			SphP[i].Density = 1000.0;
	}
#endif
		
	cout << endl << "Raytracer Init: [" << (float)timer.Time() << "] seconds, now rendering [" 
	     << Config.numFrames << "] frames:" << endl;		
		
	for(int i = 0; i < Config.numFrames; i++)
	{
		// set all frame properties based on current frame and requested keyframes
		fm->Advance();
		
		// setup camera
		Transform world2camera = fm->SetCamera();

		// recreate per frame
		VolumeIntegrator *vi = CreateVoronoiVolumeIntegrator();	
		Filter *filter       = CreateBoxFilter();
		Film *film           = CreateFilm(filter);
		Camera *camera       = NULL;
		
		if (Config.cameraFOV)	camera = CreatePerspectiveCamera(Inverse(world2camera), film);
		else                  camera = CreateOrthographicCamera(Inverse(world2camera), film);
			
		Sampler *sampler     = CreateStratifiedSampler(film, camera);
		Renderer *re         = new Renderer(sampler, camera, vi);

		// render
		if (re && scene)
			re->Render(scene,i);
			
		delete re;
	}	
	
	delete scene;
}

void arepoTestOverrides()
{
	// debugging only (Arepo2b/3b overrides)
	for (int i=0; i < NumGas; i++) {
		SphP[i].Density      = 0.1;
		SphP[i].Grad.drho[0] = 0.0;
		SphP[i].Grad.drho[1] = 0.0;
		SphP[i].Grad.drho[2] = 0.0;
	}
	
	// Arepo3b zero outer boundary
	if ( Config.filename == "test/Arepo3b" ) {
	  for(int i=0; i < NumGas; i++) {
		  if( (P[i].Pos[0] < 0.2 || P[i].Pos[0] > 0.8) ||
			    (P[i].Pos[1] < 0.2 || P[i].Pos[1] > 0.8) ||
					(P[i].Pos[2] < 0.2 || P[i].Pos[2] > 0.8) )
				SphP[i].Density = 0.0;
		}
	}
	
	// Arepo2b/3b overrides (central point)
	int dpInd;
  if ( Config.filename == "test/Arepo2b" )
		dpInd = 6;
	else if ( Config.filename == "test/Arepo3b" )
		dpInd = 45;
	else if ( Config.filename == "test/Arepo3a" )
		dpInd = 45;
	else
		terminate("probably don't want these overrides");
	
	SphP[dpInd].Density      = 20.0;
	SphP[dpInd].Grad.drho[0] = 10.0;
	SphP[dpInd].Grad.drho[1] = 10.0;
	SphP[dpInd].Grad.drho[2] = 2.0;	
}

void rtTestRenderScene(string filename)
{
	// config - camera
	Transform world2camera = LookAt(Point(Config.cameraPosition), 
	                                Point(Config.cameraLookAt), Vector(Config.cameraUp));
	
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

	// make primary classes
	VolumeIntegrator *vi  = CreateVoronoiVolumeIntegrator();
	Renderer *re          = new Renderer(sampler, camera, vi);
	VolumeRegion *vr      = NULL;
	
	// transfer function
	const Spectrum sig_a = Spectrum::FromRGB(Config.rgbAbsorb);
	TransferFunction *tf = new TransferFunction(sig_a);

	// setup transfer functions from file
	for (unsigned int i=0; i < Config.tfSet.size(); i++)
		tf->AddParseString(Config.tfSet[i]);
	
	arepoTestOverrides(); // before mesh to set Density for DTFE grads

	// voronoi mesh
	ArepoMesh *arepoMesh  = new ArepoMesh(tf);
	Scene *scene          = new Scene(vr, arepoMesh);
	
	arepoTestOverrides(); // after mesh to override Voronoi gradients for cellgrad
	
#ifdef DUMP_VORONOI_MESH
	// dump mesh in VORONOI_MESHOUTPUT format
	arepoMesh->OutputMesh();
	return;
#endif

#ifdef DUMP_MESH_TEXT
  // dump all mesh related structures in text format
	arepoMesh->DumpMesh();
	return;
#endif
	
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
		re->Render(scene,0);
			
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
	if (argc != 2 && argc != 3) {
			cout << "usage: ArepoRT <configfile.txt>\n";
			return 0;
	} else {
			// read config
			Config.ReadFile( argv[1] );

			// second argument? input snapshot and output frame number
			if( argc == 3 )
			{
			string num = argv[2];
			const char padChar = '0';

			if( num.size() < 3 ) // snapshot file number
			num.insert(0, 3-num.size(), padChar);
			size_t start_pos = Config.filename.find("NUM");
			Config.filename.replace(start_pos, 3, num);

      if( num.size() < 4 ) // output image numbering
      num.insert(0, 4-num.size(), padChar);
			start_pos = Config.imageFile.find("NUMM");
			Config.imageFile.replace(start_pos, 4, num);
			}

      if (Config.verbose)
        Config.print();

	} 
	
	// init
#ifdef ENABLE_AREPO
    Timer timer;
    timer.Start();

    Arepo arepo = Arepo(Config.filename, Config.paramFilename);

    arepo.Init(&argc,&argv);
    arepo.LoadSnapshot();

    cout << endl << "Arepo Load, Init, Meshing: [" << (float)timer.Time() << "] seconds." << endl;
#endif

	if ( Config.filename.substr(0,10) == "test/Arepo" )
		rtTestRenderScene(Config.filename);
	else
		rtIsoDiskCosmoCutoutRender();
	
	// cleanup
#ifdef ENABLE_AREPO
	arepo.Cleanup();
#endif

	return 0;
}
