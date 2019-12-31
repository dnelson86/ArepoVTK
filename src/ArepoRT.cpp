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
	
  // debugging: dump mesh in VORONOI_MESHOUTPUT format
  if(Config.dumpMeshBinary)
  {
	arepoMesh->OutputMesh();
	return;
  }

  // debugging: dump all mesh related structures in text format
  if(Config.dumpMeshText)
  {
	arepoMesh->DumpMesh();
	return;
  }

  // debugging: write SphP into stdout  
  if(Config.dumpMeshCells)
  {
	for (int i=0; i < NumGas; i++) {
			cout << "SphP[" << setw(2) << i << "] pos = " << P[i].Pos[0] << " " << P[i].Pos[1] << 
			" " << P[i].Pos[2] << " dens = " << setw(10) << SphP[i].Density 
			 << " grad.x = " << setw(10) << SphP[i].Grad.drho[0] << " grad.y = " 
			 << setw(10) << SphP[i].Grad.drho[1] << " grad.z = " << setw(10) << SphP[i].Grad.drho[2] << endl;
	}
	
	if(arepoMesh) arepoMesh->WorldBound().print("ArepoMesh WorldBound ");
  }

	// loop for requested frames	
	for(int i = Config.startFrame; i < Config.startFrame+Config.numFrames; i++)
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
			
		//vi = CreateQuadIntersectionIntegrator();

		Filter *filter       = CreateBoxFilter();
		Film *film           = CreateFilm(filter);
		Camera *camera       = CreateCamera(Inverse(world2camera), film);
			
		Sampler *sampler     = CreateStratifiedSampler(film, camera);
		Renderer *re         = new Renderer(sampler, camera, vi);

		// render
		if (re && scene)
			re->Render(scene,i);
			
		delete re;
	}	
	
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
			 << ". Author: Dylan Nelson (dnelson@mpa-garching.mpg.de)" << endl << endl;
	IF_DEBUG(cout << "DEBUGGING ENALBED.\n\n");
	
	// command line arguments
	int c;
	int curJobNum = -1;
	int expJobNum = -1;
	int snapNum = -1;
	string snapNumStr = "";
	
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
				snapNumStr = optarg;
				istringstream( snapNumStr ) >> snapNum;

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
	if( snapNumStr != "" )
	{
		const char padChar = '0';

		 // snapshot file number
		if( snapNumStr.size() < 3 )
			snapNumStr.insert(0, 3-snapNumStr.size(), padChar);
			
		size_t start_pos = Config.filename.find("NUMM");
		if( start_pos == string::npos ) {
			cout << "ERROR: Snapshot specified but no wildcard 'NUMM' in path string." << endl;
			return 0;
		}
		Config.filename.replace(start_pos, 4, snapNumStr);
		
		// second search (filename instead of directory name)
		start_pos = Config.filename.find("NUMM");
		if( start_pos != string::npos )
			Config.filename.replace(start_pos, 4, snapNumStr);

		 // output image numbering
		if( snapNumStr.size() < 4 )
			snapNumStr.insert(0, 4-snapNumStr.size(), padChar);
			
		// replace in output filename
		start_pos = Config.imageFile.find("NUMM");
		if( start_pos != string::npos ) {
			Config.imageFile.replace(start_pos, 4, snapNumStr);
		}
		else {
			cout << "ERROR: Snapshot specified but no wildcard 'NUMM' in imagePath string." << endl;
			return 0;
		}
		
		// subbox0 render: set frame number equal to snapshot number
		cout << "Subbox0 Render: setting frame number equal to snap number!" << endl;
		Config.startFrame = snapNum;
		Config.numFrames = 1;
		// end subbox0 render
	}
	
	if (Config.verbose)
		Config.print();
	
	// init
  Timer timer;
  timer.Start();

  Arepo arepo = Arepo(Config.filename, Config.paramFilename);
  arepo.Init(&argc,&argv);
  arepo.LoadSnapshot();
  arepo.ComputeQuantityBounds();

  cout << endl << "Arepo Load, Init, Meshing: [" << (float)timer.Time() << "] seconds." << endl;

	// TODO: SIGTERM
	// setup handler which recieves the SIGTERM, then sets sigTermGlobalVarFlag=1
	
	// render
	rtRenderFrames();
	
	// cleanup
	arepo.Cleanup();

	return 0;
}
