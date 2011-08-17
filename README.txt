-------------------------------------------------
--- ArepoVTK: The Arepo Visualization ToolKit ---
-------------------------------------------------
Primary Author: Dylan Nelson

Design:

-------
ArepoRT
-------
Goal:

Approach:

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
	- camera path splining in space (AnimatedTransform?)
	 
	+v0.4
  - interactive component (OpenGL)
	 - alternative/quick rendering modes
	 - navigation
	 - movie setup
	 - (single node only)
	 
 +v0.5
  - parallel (MPI)
	- large snapshot support and optimizations
	 - exchange_rays()?
  - GUI on node=0 as client	 
	 
 +v0.6
 - movie pipeline
  - frame metadata (XML/MKV container?)
 - time navigation (multiple snapshots)
 
 +v0.7
  - GPU acceleration?
	- intersection acceleration? (BVH / kdTree)
	- memory optimizations?
	- speed optimizations?
	
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
 - c++? python? (cross platform?)
 - Q: connect to a "server" or spawn the GUI of ArepoRT at MPI node==0 

----------------------
Data Format Conversion
----------------------
 - ParaView / VisIT
 - VTK
 - tetra, block amr, or unigrid cubes
 