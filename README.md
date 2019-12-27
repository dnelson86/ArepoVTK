
ArepoVTK: The Arepo Visualization ToolKit
=========================================

Produce high quality, presentation-ready visualizations of hydrodynamic simulations run with the [AREPO](https://arepo-code.org/) code. Ray cast through linearly reconstructed scalar and vector fields defined on an unstructured Voronoi tessellation of space. Support multi-dimensional transfer functions to investigate fluid quantities, and explore novel visualization techniques for combining such a volume rendering approach with coincident point particle sets (both luminous and dark).

A high level renderering framework is implemented in C++, and we dynamically link with AREPO to use its existing functionality to load snapshots, initialize fluid and particle data, and construct the Voronoi mesh and its connectivity structures. Transfer function design is the burden of the user and assumes an expert knowledge of the data present in the snapshots. This is specified, along with all other rendering options, in a configuration file read at run time. 

### Installation

First, make sure you have a recent C++ compiler (gcc or intel). On a standard cluster, load a set of required modules by executing `module load gcc gsl fftw-serial hdf5-serial impi`. If you are on a laptop or otherwise don't have the `module` command available, you must install these libraries (consult the [AREPO user documentation](https://gitlab.mpcdf.mpg.de/vrs/arepo/wikis/userguide/Getting%20started) for more details). 

Next, download the ArepoVTK source as well as the public AREPO source:

```bash
git clone https://github.com/dnelson86/ArepoVTK.git
cd ArepoVTK/
git clone https://github.com/dnelson86/arepo
```

Install `libpng` if it is not already available on your system:

```bash
git clone git://git.code.sf.net/p/libpng/code libpng
cd libpng
./configure
make
```

Build the `libarepo.a` shared library used within ArepoVTK:

```bash
cd ../arepo
make libarepo.a -j
```

Finally, build ArepoVTK itself:

```bash
cd ..
make -j
```

Executing `./ArepoRT` with no options should produce the following output:

```
      ___           ___           ___           ___           ___                    ___           ___
     /\  \         /\  \         /\  \         /\  \         /\  \                  /\  \         /\  \
    /::\  \       /::\  \       /::\  \       /::\  \       /::\  \                /::\  \        \:\  \
   /:/\:\  \     /:/\:\  \     /:/\:\  \     /:/\:\  \     /:/\:\  \              /:/\:\  \        \:\  \
  /::\~\:\  \   /::\~\:\  \   /::\~\:\  \   /::\~\:\  \   /:/  \:\  \            /::\~\:\  \       /::\  \
 /:/\:\ \:\__\ /:/\:\ \:\__\ /:/\:\ \:\__\ /:/\:\ \:\__\ /:/__/ \:\__\          /:/\:\ \:\__\     /:/\:\__\
 \/__\:\/:/  / \/_|::\/:/  / \:\~\:\ \/__/ \/__\:\/:/  / \:\  \ /:/  /          \/_|::\/:/  /    /:/  \/__/
      \::/  /     |:|::/  /   \:\ \:\__\        \::/  /   \:\  /:/  /              |:|::/  /    /:/  /
      /:/  /      |:|\/__/     \:\ \/__/         \/__/     \:\/:/  /               |:|\/__/     \/__/
     /:/  /       |:|  |        \:\__\                      \::/  /                |:|  |
     \/__/         \|__|         \/__/                       \/__/                  \|__|

   v0.44 (Dec 25 2019). Author: Dylan Nelson (dnelson@mpa-garching.mpg.de)

Usage: ArepoRT <configfile.txt> [-s snapNum] [-j jobNum] [-e expandedJobNum]

```


### Getting Started

After compilation succeeds, run the first basic test:

```bash
./ArepoRT tests/test2_config.txt
```

this should produce a 600x600 pixel image `test_frame2.png` which is identical to the following:

![ArepoVTK test2 result](tests/test_frame2.png?raw=true "ArepoVTK test2 result")

This is a result of a simple ray-tracing through a small test setup of a uniform grid of 2^3 (i.e. eight) cells within a box [0,1], each with uniform mass. A ninth, central point at [0.5, 0.5, 0.5] is inserted with higher mass. As a result, the gas density field peaks in the center and falls off radially, modulo the imprint of the tessellation geometry.

The transfer function is defined in the configuration file: in this case, there is only one TF added and it is given by the string `constant_table Density idl_33_blue-red 0.5 20`. The 'constant' means that we will add light to a ray in linear proportion to the 'Density' field it samples. The color of the light is sampled from a 'table', which is specified by its name 'idl_33_blue-red'. This colormap is stretched between a minimum Density value of '0.5' and a maximum of '20' (code units), and outside of this range gas will not contribute to a ray.

Next, let's run a permutation of this test on the same "simulation snapshot":

```bash
./ArepoRT tests/test2b_config.txt
```

which should produce the image `test_frame2b.png` as shown below:

![ArepoVTK test2b result](tests/test_frame2b.png?raw=true "ArepoVTK test2b result")

Several configuration choices were changed from the first image. First, the camera was moved such that it views the simulation domain from an oblique angle, rather than directly 'face-on'. The geometry of the bounding box and the single octagonal Voronoi cell in the center of the domain is clear. Second, we have changed the transfer function to `constant Density 1.0 0.2 0.0` which is even simpler than above: a fixed color specified by the RGB triplet {R=1, G=0.2, B=0} (i.e. orange) is added to a ray each time it samples gas 'Density', weighted by the value of that density. Finally, we have changed the `viStepSize` parameter from zero to `0.05`. If `viStepSize = 0`, then ArepoVTK samples each Voronoi cell exactly once, at the midpoint of the line segment defined by the two intersection points of a ray as it enters and exits that cell. On the other hand, if `viStepSize > 0` as in the second example, we take strictly step along each ray and take equally spaced samples 0.05 (world space, i.e. code units) apart.

Third, the test:

```bash
./ArepoRT tests/test3_config.txt
```

explores a slighly larger case. This is a 4^3 uniform grid of cell generating sites, again within a simulation volume defined as [0,1], with a single centered cell added, for 65 total cells.


### Gallery of Examples

All of the gas images of the [Illustris Explorer](https://www.illustris-project.org/explorer/) were generated with ArepoVTK, using the natural neighbor interpolation (NNI) method. The configuration files are available under `tests/illustris_box*`.

* [The Universe in Gas](https://vimeo.com/77612968) video was made with ArepoVTK, showing gas iso-density and iso-temperature contours within a 20 Mpc/h cosmological volume, each using a set of `gaussian_table` transfer functions. The configuration files are available under `tests/cosmoRot*`.

* [Stirring Coffee with a Spoon in 3D](https://vimeo.com/72435369) video was made with ArepoVTK, using gaussian transfer functions on gas density. It is made up of three sequences: the initial rotation, the time-evolving sequence, and the ending. The configuration files are available under `tests/spoon*`.

Figure 13 of [Nelson et al. (2016)](https://arxiv.org/abs/1503.02665) shows gas iso-temperature contours around a single galaxy halo. The configuration files are available under `tests/zoom_Nelson16*`.

A gallery of test images during development is available on the [ArepoVTK Development Gallery](https://wwwmpa.mpa-garching.mpg.de/~dnelson/#arepovtkgal).


### Acknowledgments, Contributing, and Contact

The principal author of ArepoVTK Dylan Nelson (dnelson@mpa-garching.mpg.de). If you use this code in a publication, please cite [Nelson et al. (2013)](http://ui.adsabs.harvard.edu/abs/2013MNRAS.429.3353N).

Please contact us with any questions or comments. If you make any changes, updates, fixes, or improvements to this code, please consider making your efforts available to the broader community by contributing back to this project. The most effective way to do so is by submitting an Issue or Pull Request on github.


### Version Roadmap

Current Version: 0.44.alpha1

 +v0.1
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
