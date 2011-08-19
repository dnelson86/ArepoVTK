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

void rtTestRenderScene(string filename)
{
		// config - camera
		Point cameraPos    = Point(10.5,10.5,0.0); //centered on z=0 plane edge, x/y axis aligned
		Point cameraLook   = Point(10.5,10.5,1.0);
		
		//Point cameraPos      = Point(10.0,10.5,40.0); //centered on x=1.0 plane edge
		//Point cameraLook     = Point(10.5,10.5,40.5); //looking -zhat, y/z axis aligned
		
		//Point cameraPos    = Point(12.0,12.0,38.0); //angled above
		//Point cameraLook   = Point(10.5,10.5,40.5); //looking at center of box
		
		Vector cameraVecUp = Vector(0.0,1.0,0.0);
		
		// set transform
		Transform volume2world;
		Transform world2camera;

		volume2world = volume2world * Translate(Vector(10.0,10.0,40.0));
		world2camera = world2camera * LookAt(cameraPos, cameraLook, cameraVecUp);
		
		IF_DEBUG(world2camera.print("world2camera:"));
		IF_DEBUG(Inverse(world2camera).print("camera2world:"));
		IF_DEBUG(volume2world.print("volume2world:"));
		
		// Camera (Filter, Film)
		Filter *filter       = CreateBoxFilter();
		Film *film           = CreateFilm(filter);
		Camera *camera       = CreateOrthographicCamera(Inverse(world2camera), film);
		
		// Sampler
		Sampler *sampler     = CreateStratifiedSampler(film, camera);
		
		// Integrator
		//VolumeIntegrator *vi = CreateEmissionVolumeIntegrator(Config.viStepSize);
		VolumeIntegrator *vi  = CreateVoronoiVolumeIntegrator();
		
		// renderer
		Renderer *re         = new Renderer(sampler, camera, vi);
		
		// transfer function
		TransferFunction *tf = CreateTransferFunction();
		
		// debugging:
		Spectrum s1 = Spectrum::FromRGB(Config.rgbEmit);
		Spectrum s2 = Spectrum::FromNamed("green");
		Spectrum s3 = Spectrum::FromNamed("blue");
		tf->AddConstant(TF_VAL_DENS,s2);
		//tf->AddTophat(TF_VAL_DENS,5.0,10.0,spec);
		//tf->AddGaussian(TF_VAL_DENS,0.0,15.0,4.0,0.2,s1);
		//tf->AddGaussian(TF_VAL_DENS,0.0,15.0,8.0,0.2,s2);
		
		// create volume/density/scene geometry (debugging only)
		//VolumeRegion *vr     = CreateGridVolumeRegion(volume2world, filename);
		VolumeRegion *vr      = NULL;
		
		// voronoi mesh
		ArepoMesh *arepoMesh  = CreateArepoMesh(tf, volume2world);

		// debugging only (Arepo2b overrides)
		for (int i=0; i < N_gas; i++) {
				SphP[i].Density      = 0.05;
				SphP[i].Grad.drho[0] = 0.0;
				SphP[i].Grad.drho[1] = 0.0;
				SphP[i].Grad.drho[2] = 0.0;
				if (i == 6) {
						SphP[i].Density      = 10.0;
						SphP[i].Grad.drho[0] = 0.2;
						SphP[i].Grad.drho[1] = 0.2;
						SphP[i].Grad.drho[2] = 0.0;
				}
				cout << "SphP[" << setw(2) << i << "] dens = " << setw(10) << SphP[i].Density 
						 << " grad.x = " << setw(10) << SphP[i].Grad.drho[0] << " grad.y = " 
						 << setw(10) << SphP[i].Grad.drho[1] << " grad.z = " << setw(10) << SphP[i].Grad.drho[2] << endl;
	  }
		
		// scene
		Scene *scene          = new Scene(vr, arepoMesh);
		
#ifdef DEBUG
    if(arepoMesh) arepoMesh->WorldBound().print("ArepoMesh WorldBound ");
		if(vr)        vr->WorldBound().print("VolumeRegion WorldBound ");
		
		//Point testp(0.2, 0.2, 0.2);
		//cout << "VR Test @ (" << testp.x << " " << testp.y << " " << testp.z << ") result = "
		//		 << vr->Density(testp) << endl;
		
		float vals[2];
		vals[1] = 0.0; //utherm
		
		for (float rho = 0.0; rho < 10.0; rho += 1.0) {
		  vals[0] = rho;
			Spectrum Lve = tf->Lve(vals);
			
			cout << "TF Test @ (" << vals[1] << "," << vals[0] << ") R = " << setw(7) << Lve.r() 
					 << " G = " << setw(7) << Lve.g() << " B = " << setw(7) << Lve.b() << endl;
		}
#endif
				 
		// render
		if (filter && film && camera && sampler && vi && scene && tf && re)
				re->Render(scene);
				
		delete re;
		delete scene;
}

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
