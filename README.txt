-------------------------------------------------
--- ArepoVTK: The Arepo Visualization ToolKit ---
-------------------------------------------------

Primary Author: Dylan Nelson (dnelson@cfa.harvard.edu)

Current Version: 0.3 (22 Aug 2011)

Installation/Compilation:

 (1) Get ArepoVTK Source:   hg clone /n/hernquistfs1/hgrepos/ArepoVTK/ ~/ArepoVTK/
 (2) Get Arepo Source:      svn checkout http://www.mpa-garching.mpg.de/svn/cosmo-group/Arepo ~/Arepo/
 (3) Copy Arepo Config:     cp ~/ArepoVTK/test/Config_ArepoVTK.sh ~/Arepo/
 (4) Build Arepo:           cd ~/Arepo && make CONFIG=Config_ArepoVTK.sh EXEC=Arepo
 (5) Build ArepoVTK:        cd ~/ArepoVTK && make clean && make
 (6) Run ArepoRT Test:      ./ArepoRT

Design:

-------
ArepoRT
-------
Goal:     Produce high quality, presentation-ready visualizations of hydrodynamic simulations run 
          with Arepo. Ray cast through linearly reconstructed scalar and vector fields defined on an 
			    unstructured Voronoi tessellation of space. Support multi-dimensional transfer functions to 
			    investigate fluid quantities, and explore novel visualization techniques for combining such 
			    a volume rendering approach with coincident point particle sets (both luminous and dark).

Approach: Todo.

Version Roadmap:
 + v0.1
  - program shell and commandline functionality
  - radiative transfer framework:
	 - pts, vector, matrices, rays, bboxes
	 - camera: ortho projection, LookAt
	 - stratified (no jitter or MC) sampler
	 - volume emission (no scattering) integrator
	 - renderer, task/image plane subdivision
	- image output: raw text floats / TGA
	- volumegrid input: "debug" scene descriptor format
	
 +v0.2
  - interface with Arepo or cut needed portions to acquire tesselation
	- wireframe T. rendering
	
 +v0.3
  - voronoi ray casting
	 - cell averaged constant
	 - gradients
	- transfer function on a single quantity
  - configuration file
	
+v0.4
  - parallel (MPI)
	- large snapshot support and optimizations
   - use treefind for entry cells
   - handle local and foreign ghosts
	 - exchange_rays() type approach
	 
+v0.6
 - camera path splining in space (AnimatedTransform?)
 - movie pipeline
  - frame metadata (XML/MKV container?)
 - time navigation (multiple snapshots, interpolation?)	 
	 
+v0.5
  - interactive component (OpenGL)
	 - alternative/quick rendering modes
	 - navigation
	 - movie setup
	 - (single node only)
  - GUI on node=0 as client (openWindow=true)

 +v0.7
  - GPU acceleration (ray tracing)?
	- intersection acceleration? (BVH / kdTree)
	- memory optimizations (wipe out / rearrange some Arepo stuff)
	- speed optimizations?
	 - scratch Volume2World transform?
	
 +v0.8
  - tetrahedra decomposition: VisIT exporter
	- light sources (star particles)
	- single scattering volume integrator (MC)
	- DOF/motion blur

--------------------------
Arepo Interactive/Explorer
--------------------------
Goal:

Ideas:
 - lightweight desktop client
 - realtime and interactive
 - openGL
  - windowing: freeglut/glew, QT?
 - c++? python? (cross platform?)
 - Q: connect to a "server" or spawn the GUI of ArepoRT at MPI node==0 

----------------------
Data Format Conversion
----------------------
 - ParaView / VisIT
 - VTK
 - tetra, block amr, or unigrid cubes
 
