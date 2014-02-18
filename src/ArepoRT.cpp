/*
 * ArepoRT.cpp
 * dnelson
 */
 
#include <unistd.h>
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

void rtRenderFrames()
{
	Timer timer;
	timer.Start();
	
	const Spectrum sig_a = Spectrum::FromRGB(Config.rgbAbsorb);

	// constant setup		
	TransferFunction *tf  = new TransferFunction(sig_a);
	FrameManager *fm      = new FrameManager(Config.kfSet);
	
	// sampling based on tree only, or full mesh construction?
	ArepoTree *arepoTree = NULL;
	ArepoMesh *arepoMesh = NULL;
	
	if( Config.nTreeNGB )
		arepoTree = new ArepoTree(tf);
	else
		arepoMesh = new ArepoMesh(tf);
	
	Scene *scene = new Scene(NULL, arepoMesh, arepoTree);
	
	// setup transfer function
	for (unsigned int i=0; i < Config.tfSet.size(); i++)
		tf->AddParseString(Config.tfSet[i]);

	cout << endl << "Raytracer Init: [" << (float)timer.Time() << "] seconds, now rendering [" 
	     << Config.numFrames << "] frames:" << endl;		
	
	// loop for requested frames	
	for(int i = Config.startFrame; i < Config.numFrames; i++)
	{
		// set all frame properties based on current frame and requested keyframes
		fm->Advance(i);
	
		// setup camera
		Transform world2camera = fm->SetCamera();
		VolumeIntegrator *vi = NULL;
		
		// recreate per frame
		if( Config.nTreeNGB )
			vi = CreateTreeSearchVolumeIntegrator();
		else
			vi = CreateVoronoiVolumeIntegrator();

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
	
	for( int i=0; i < NumGas; i++ )
		if( SphP[i].OldMass < 0.0 ) {
			cout << "Error: After finished, found [" << i << "] entropy = " << SphP[i].OldMass << endl;
			exit(1207);
		}
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
	Scene *scene          = new Scene(vr, arepoMesh, NULL);
	
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

	cout << "      ___           ___           ___           ___           ___                    ___           ___     \n"
       << "     /\\  \\         /\\  \\         /\\  \\         /\\  \\         /\\  \\                  /\\  \\         /\\  \\    \n"
			 << "    /::\\  \\       /::\\  \\       /::\\  \\       /::\\  \\       /::\\  \\                /::\\  \\        \\:\\  \\   \n"
			 << "   /:/\\:\\  \\     /:/\\:\\  \\     /:/\\:\\  \\     /:/\\:\\  \\     /:/\\:\\  \\              /:/\\:\\  \\        \\:\\  \\  \n"
			 << "  /::\\~\\:\\  \\   /::\\~\\:\\  \\   /::\\~\\:\\  \\   /::\\~\\:\\  \\   /:/  \\:\\  \\            /::\\~\\:\\  \\       /::\\  \\ \n"
			 << " /:/\\:\\ \\:\\__\\ /:/\\:\\ \\:\\__\\ /:/\\:\\ \\:\\__\\ /:/\\:\\ \\:\\__\\ /:/__/ \\:\\__\\          /:/\\:\\ \\:\\__\\     /:/\\:\\__\\\n"
			 << " \\/__\\:\\/:/  / \\/_|::\\/:/  / \\:\\~\\:\\ \\/__/ \\/__\\:\\/:/  / \\:\\  \\ /:/  /          \\/_|::\\/:/  /    /:/  \\/__/\n"
			 << "      \\::/  /     |:|::/  /   \\:\\ \\:\\__\\        \\::/  /   \\:\\  /:/  /              |:|::/  /    /:/  /     \n"
			 << "      /:/  /      |:|\\/__/     \\:\\ \\/__/         \\/__/     \\:\\/:/  /               |:|\\/__/     \\/__/      \n"
			 << "     /:/  /       |:|  |        \\:\\__\\                      \\::/  /                |:|  |                  \n"
			 << "     \\/__/         \\|__|         \\/__/                       \\/__/                  \\|__|                  ";


	cout     << endl << endl
	     << "   v" << AREPO_RT_VERSION << " (" << __DATE__      << ")"
			 << ". Author: Dylan Nelson (dnelson@cfa.harvard.edu)" << endl << endl;
	IF_DEBUG(cout << "DEBUGGING ENALBED.\n\n");
	
	// command line arguments
	int c;
	int curJobNum = -1;
	int expJobNum = -1;
	string snapNum = "";
	
	opterr = 0;
	
	while ((c = getopt( argc, argv, "s:j:e:")) != -1 )
		switch (c)
		{
			case 'j':
			{
				string num = optarg;
				istringstream( num ) >> curJobNum;
				
				break;
			}
			case 'e':
			{
				string num = optarg;
				istringstream( num ) >> expJobNum;

				break;
			}
			case 's':
			{
				string snapNum = optarg;

				break;
			}
			default:
				cout << "When does this happen?" << endl;
				return 0;
		}
	
	// should be one non-option, the configuration file parameter
	if( optind != argc - 1 ) {
		cout << "Usage: ArepoRT <configfile.txt> [-s snapNum] [-j jobNum] [-e expandedJobNum]" << endl << endl;
		return 0;
	}
	
	// read config
	Config.ReadFile( argv[optind] );
	
	// process 'j' input
	if( curJobNum != -1 )
	{
		Config.curJobNum = curJobNum;
		// require that the current job number is reasonable given the total number of jobs
		if( curJobNum < 0 || curJobNum >= Config.totNumJobs ) {
			cout << "ERROR: Specified curJobNum [" << curJobNum << "] but Config.totNumJobs is [" << Config.totNumJobs << "]" << endl;
			return 0;
		} else {
			cout << "Executing job number [" << curJobNum+1 << "] of [" << Config.totNumJobs << "]." << endl;
		}
	} else {
		// require that we have not specified the current job number in this case
		if( Config.totNumJobs > 0 ) {
			cout << "Error: Have totNumJobs > 0 but the current job number is unspecified!" << endl;
			return 0;
		}
		if( expJobNum > 0 ) {
			cout << "Error: Have expandedJobNum > 0 but the current job number is unspecified!" << endl;
			return 0;
		}
	}

	// process 'e' input
	if( expJobNum != -1 )
	{
		Config.expandedJobNum = expJobNum;
		// require that the expanded job number is reasonable
		if( expJobNum < 0 || expJobNum >= Config.totNumJobs*pow(Config.jobExpansionFac,2) ) {
			cout << "ERROR: Specified expJobNum is out of bounds." << endl;
			return 0;
		} else {
			cout << " -- -- Executing expanded job number [" << expJobNum+1 << "] of ["
			     << Config.totNumJobs*pow(Config.jobExpansionFac,2) << "]." << endl;
		}
	} else {
		// require that we have not specified a jobExpansionFactor above one in this case
		if( Config.jobExpansionFac != 1 ) {
			cout << "Error: Have jobExpansionFac above unity, but expandedJobNumber not specified." << endl;
			return 0;
		}
	}
	
	// process 's' input
	if( snapNum != "" )
	{
		const char padChar = '0';

		 // snapshot file number
		if( snapNum.size() < 3 )
			snapNum.insert(0, 3-snapNum.size(), padChar);
		size_t start_pos = Config.filename.find("NUM");
		if( start_pos == string::npos ) {
			cout << "ERROR: Snapshot specified but no wildcard 'NUM' in path string." << endl;
			return 0;
		}
		Config.filename.replace(start_pos, 3, snapNum);

		 // output image numbering
		if( snapNum.size() < 4 )
			snapNum.insert(0, 4-snapNum.size(), padChar);
			
		Config.imageFile += "_" + snapNum;
	}
	
	if (Config.verbose)
		Config.print();
	
	// init
  Timer timer;
  timer.Start();

  Arepo arepo = Arepo(Config.filename, Config.paramFilename);
  arepo.Init(&argc,&argv);
  arepo.LoadSnapshot();

  cout << endl << "Arepo Load, Init, Meshing: [" << (float)timer.Time() << "] seconds." << endl;

	// TODO: SIGTERM
	// setup handler which recieves the SIGTERM, then sets sigTermGlobalVarFlag=1
	
	// render
	if ( Config.filename.substr(0,10) == "test/Arepo" )
		rtTestRenderScene(Config.filename);
	else
		rtRenderFrames();
	
	// cleanup
	arepo.Cleanup();

	return 0;
}
