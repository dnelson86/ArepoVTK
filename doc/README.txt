-------------------------------------------------
--- ArepoVTK: The Arepo Visualization ToolKit ---
-------------------------------------------------

Primary Author: Dylan Nelson (dnelson@cfa.harvard.edu)

Current Version: 0.43 (18 Feb 2014)

Installation/Compilation:

 (1) Get ArepoVTK Source:   hg clone /n/hernquistfs1/hgrepos/ArepoVTK/ ~/ArepoVTK/
 (2) Get Arepo Source:      svn checkout http://www.mpa-garching.mpg.de/svn/cosmo-group/Arepo ~/ArepoVTK/Arepo/
 (3) Copy Config:           cp ~/ArepoVTK/Config_ArepoVTK.sh ~/ArepoVTK/Arepo/Config.sh
 (3) Build libarepo.a:      cd ~/ArepoVTK/Arepo && make libarepo.a
 (4) Build ArepoVTK:        cd ~/ArepoVTK && make clean && make
 (5) Run ArepoRT Test:      ./ArepoRT test/config.txt

Notes:

 * To build libarepo.a, delete main() from main.c and compile Arepo (ignore
 * missing main error), then build libarepo.a.

 * All structure is copied into local allvars.h and must match compile options
 * and current revision in Arepo allvars.h, otherwise All.BoxSize likely zero
 * and rays likely all terminated at their start.

-------
ArepoRT
-------
Goal:     Produce high quality, presentation-ready visualizations of hydrodynamic simulations run 
          with Arepo. Ray cast through linearly reconstructed scalar and vector fields defined on an 
          unstructured Voronoi tessellation of space. Support multi-dimensional transfer functions to 
          investigate fluid quantities, and explore novel visualization techniques for combining such 
          a volume rendering approach with coincident point particle sets (both luminous and dark).

Approach: A high level renderering framework is implemented in c++ following the design methodology 
          in the book "Physically Based Rendering" (2ed) by M. Pharr and G. Humphreys. Integration 
          with Arepo follows the strategy of SUNRISE (written by P. Jonsson) - dynamically linking 
          with Arepo (written by V. Springel) and using its
          existing functionality to load a snapshot, initialize fluid and particle data, and construct
          the mesh and its connectivity structures. We borrow heavily from voronoi_makeimage_new.c to
          raytrace through the Voronoi mesh. Transfer function design is the burden of the user and
          assumes an expert knowledge of the data present in the snapshots. These are specified, 
          along with all other rendering options, in a configuration file read at run time. 

Usage:    See the included user manual for details.

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
  - interface with Arepo
   - init and snapshot loading
   - acquire mesh construction and connectivity
    - wireframe tetra/bounding box rendering
   - wu's line algorithm
	
 +v0.3
  - voronoi ray casting
   - cell averaged constant
   - gradients
  - attenuation via optical depth along LOS
   - transfer function on a single quantity
  - configuration file

 +v0.35
  - large snapshot support and optimization
   - use treefind for entry cells
   - handle local ghosts
   - 128^3 cosmo rendering
  - parallel (threads) on shared memory node

 +v0.39
  - external colortables
  - voronoi cell algorithm
  - SPH kernel and IDW of natural neighbors interp methods
  - tetra mesh walking
   - DTFE 
  - use internal DC connectivity instead of alternative Sunrise

 +v0.4
  - multiple frames on single snapshot
  - camera path splining/tweening in space, orbits
  - fix NNI auxMesh allocation per task instead of per thread

 +v0.41
  - support for SphKernel and IDW using N-nearest neighbor spheres
    using NGBTree searches, skip all mesh related tasks including
    construction and memory overhead

 +v0.42
  - calculate raw line integrals (do not apply TF) in addition to images
  - goal: rendering 1820^3 box
   - custom read_ic() with selective loading of only particles in camera frustrum
   - split camera/image into subchunks, each handled by a separate instance
   - maskfile approach to precompute which particles each instance needs to load
 
 +v0.43
  - expanded job subdivision
  - adaptive spatial stepping with tree based ray tracing
 
 +v0.44
  - within custom loading, load DM, Stars, and BHs for rendering (one at a time)
   - transfer functions fields for each that make sense for Illustris
  - fisheye (180 deg) and environment (360 deg) non-linear cameras
  - QuadIntersectionIntegrator (based on single-pass wireframe rendering)
  - remove RasterizeLine approach, do bbox lines inside rays based on bbox intersection
   - attenuation based on ray optical depth before reaching line

 +v0.45
  - visualization voronoi mesh based within rays upon each edge intersection
  - "splat" integrator via disc intersection testing + kdtree?
  - stars: use BC03 stellar photometric tables for U,B,V,K,g,r,i,z band images

 +v0.46
  - keyframe transfer function settings
  - time navigation (multiple snapshots, time interpolation)
   - subboxes: render bbox's, update from disk at different interval than fullbox
  - movie pipeline
   - frame metadata (XML/MKV container?)  

 +v0.47
  - fix/verify exact NNI
  - Watson-Sambridge NNI
  - Liang-Hale NNI

 +v0.48
  - load group catalogs, merger trees
   - track halos in time

 +v0.49
  - 2D transfer functions, e.g. f(rho,T)	 
  - derived fields for TFs (e.g. temp, coolingrate)
  - johnston convolved planck BB temp emission TF

 +v0.50
  - custom memory downsizing (minimize P,SphP)
  - robust restart functionality
  - load multiple particle types simultaneously: gas/dm, gas/stars
   - multiple meshes (e.g. allow DM NNI and gas NNI in TF) 

 +v0.55
  - parallel (MPI) on distributed memory cluster
   - handle foreign ghosts
   - exchange_rays() type approach

 +v0.6
  - interactive component (OpenGL)
   - alternative/quick rendering modes (splatting)
   - navigation
   - movie setup
   - (single node only)
  - GUI on node=0 as client (openWindow=true)

 +v0.61
  - GPU acceleration (ray tracing)?
   - intersection acceleration? (BVH / kdTree)
   - memory optimizations (wipe out / rearrange some Arepo stuff)
   - speed optimizations?
  - discrete NNI on GPU (sub-blocks to meeting sampling requirement)
	
 other
  - light sources (star particles)?
  - single scattering volume integrator (MC)?
  - DOF/motion blur?

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
 
