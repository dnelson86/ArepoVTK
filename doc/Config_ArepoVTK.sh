#!/bin/bash            # this line only there to enable syntax highlighting in this file

##################################################
#  Enable/Disable compile-time options as needed #
##################################################


#--------------------------------------- Basic operation mode of code

PERIODIC
#TWODIMS
#AXISYMMETRY                    # This is for axisymmetry in cylindrical coordinates (requires TWODIMS and a stationary mesh)
#ONEDIMS
#LONG_X=10.0
#LONG_Y=2.0
#LONG_Z=10.0
#REFLECTIVE_X=1 #=2            # if set to 2, the boundary is inflow/outflow */
#REFLECTIVE_Y=1 #=2
#REFLECTIVE_Z=1 #=2

#MHD
#MHD_DIVBCLEANING
#MHD_POWELL
#MHD_SEEDFIELD
#COOLING
#USE_SFR
#SFR_KEEP_CELLS
#SINKS
#GAMMA=1.4
#ISOTHERM_EQS
#USE_ENTROPY_FOR_COLD_FLOWS
#ENTROPY_MACH_THRESHOLD=1.1
#PREHEATING

#----------------------------------------MPI/Threading Hybrid
#NUM_THREADS=4
#THREAD_SAFE_COSTS                        #exact cost for domain decomposition, but slow in tree walk and in domain decomposition cost sumup
#SPIN_LOCK                               #use spin lock instead of mutex  

#--------------------------------------- Mesh Type
#AMR
VORONOI

#--------------------------------------- Riemann solver
#VARIABLE_GAMMA
#RIEMANN_HLLC
#RIEMANN_ROSUNOV
#RIEMANN_HLLD
#RIEMANN_GAMMA

#--------------------------------------- Reconstruction
#TVD_SLOPE_LIMITER
#DISABLE_TIME_EXTRAPOLATION              /* use only when you know exactly what you are doing. activating this option will make your results wrong but can tell you about the behaviour of your code */
#DISABLE_SPATIAL_EXTRAPOLATION           /* use only when you know exactly what you are doing. activating this option will make your results wrong but can tell you about the behaviour of your code */

#--------------------------------------- Mesh motion and regularization
#VORONOI_STATIC_MESH
REGULARIZE_MESH_CM_DRIFT
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
#REGULARIZE_MESH_FACE_ANGLE
#OUTPUT_MESH_FACE_ANGLE
#CALCULATE_VERTEX_VELOCITY_DIVERGENCE

#--------------------------------------- Time integration options
#FORCE_EQUAL_TIMESTEPS    # this chooses a variable but global timestep
TREE_BASED_TIMESTEPS


#--------------------------------------- Direct Volume Raycast Rendering 
#DVR_RENDER=0                             # direct volumetric raycasting, 0->stand-alone, 1->on-the-fly
#DVR_RENDER_SMOOTH                        # smooth field 
#DVR_RENDER_ORTHOGONAL                    # orthogonal projection 


#--------------------------------------- Image generation
#VORONOI_MESHOUTPUT                      # 2D and 3D mesh output
#VORONOI_IMAGES_FOREACHSNAPSHOT
#VORONOI_FREQUENT_IMAGES                 # creates images with frequency 'TimeBetweenImages' given in parameterfile, independent of snapshots  
#VORONOI_FIELD_DUMP_PIXELS_X=1536
#VORONOI_FIELD_DUMP_PIXELS_Y=150
#VORONOI_VELOCITY_FIELD_2D
#VORONOI_FIELD_COMPENSATE_VX=4.0
#VORONOI_FIELD_COMPENSATE_VY=0
#VORONOI_NEW_IMAGE
#VORONOI_PROJ_TEMP                            #project T instead of u
#VORONOI_PROJ                             # do projection along any predefined direction
#VORONOI_MULTIPLE_PROJECTIONS             # do face-on and edge-on projections by swapping y and z axes



#--------------------------------------- Refinement and derefinement
#REFINEMENT_SPLIT_CELLS
#REFINEMENT_MERGE_CELLS
#REFINEMENT_MERGE_PAIRS
#REFINEMENT_VOLUME_LIMIT
#REFINEMENT_HIGH_RES_GAS
#DEREFINE_ONLY_DENSE_GAS
#NODEREFINE_BACKGROUND_GRID
#DEREFINE_GENTLY

#--------------------------------------- Mesh-relaxing or mesh-adding (this will not carry out a simulation)
#MESHRELAX                     # this keeps the mass constant and only regularizes the mesh
MESHRELAX_DENSITY_IN_INPUT
#ADDBACKGROUNDGRID=16

#--------------------------------------- Gravity treatment
#SELFGRAVITY                   # switch on for self-gravity     
#NO_GAS_SELFGRAVITY            # switch off gas self-gravity in tree 
#GRAVITY_NOT_PERIODIC          # if gravity is not to be treated periodically
#EXACT_GRAVITY_FOR_PARTICLE_TYPE=4 #N-squared fashion gravity for a small number of particles of the given type
#EXTERNALGRAVITY               # switch on for external potential
#EXTERNALGY=0.0
#EXTERNALDISKPOTENTIAL
#FIXED_GRAVITATIONAL_SOFTENINGS_FOR_CELLS
#ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS
#ENFORCE_JEANS_STABILITY_OF_CELLS    # this imposes an adaptive floor for the temperature
#EVALPOTENTIAL                 # computes gravitational potential
#EXTERNALSHEETY
#COMPUTE_POTENTIAL_ENERGY

#--------------------------------------- TreePM Options
#PMGRID=128
#ASMTH=1.25
#RCUT=6.0

#PLACEHIGHRESREGION=2
#ENLARGEREGION=1.1
#GRIDBOOST=2
#ONLY_PM
#PMPERIODIC_LOWMEM_THREADED       # This replaces the standard pm_periodic.c with a version that uses less peak memory and is slightly faster
                                  # but this routine only works well for homogeneously sampled boxes. It can also be used with threads, in this case both OPENMP NUM_THREADS should be set


#--------------------------------------- Things that are always recommended
#AUTO_SWAP_ENDIAN_READIC        # Enables automatic ENDIAN swapping for reading ICs
#PEANOHILBERT_EXTEND_DYNAMIC_RANGE
CHUNKING                 # will calculated the gravity force in interleaved blocks. This can reduce imbalances in case multiple iterations due to insufficient buffer size need to be done
                          

#---------------------------------------- Single/Double Precision
#DOUBLEPRECISION=1
#DOUBLEPRECISION_FFTW
#OUTPUT_IN_DOUBLEPRECISION # snapshot files will be written in double precision
#INPUT_IN_DOUBLEPRECISION  # initial conditions are in double precision

#---------------------------------------- On the fly FOF groupfinder 
#FOF                                # enable FoF output
#FOF_PRIMARY_LINK_TYPES=2           # 2^type for the primary dark matter type
#FOF_SECONDARY_LINK_TYPES=1+16+32   # 2^type for the types linked to nearest primaries
#FOF_GROUP_MIN_LEN=32               # default is 32
#FOF_LINKLENGTH=0.16                # Linkinglength for FoF (default=0.2)
#FOF_FUZZ_SORT_BY_NEAREST_GROUP=0   # sort fuzz particles by nearest group and generate offset table in catalog (=1 writes nearest group number to snapshot)
#FOF_STOREIDS                       # Stores IDs in group/subfind catalogue (in principle redundant if group catalogues are produced on the fly as then the snapshots are reshuffled automatically)

#---------------------------------------- Subfind
#SUBFIND                            # enables substructure finder
#SAVE_HSML_IN_SNAPSHOT              # this will store hsml and density values in the snapshot files 

#SUBFIND_MEASURE_H2MASS             # special measuremenat option for mass in molecular hydrogen
  
#--------------------------------------- SFR/feedback model
#METALS
#MIN_METALLICITY_ON_STARTUP
#STELLARAGE

#SOFTEREQS
#SLOW_RELAX_TO_EOS

#-------------------------------------- AGN stuff
#BLACK_HOLES             # enables Black-Holes (master switch)
#BONDI                   # Bondi-Hoyle style accretion model
#ENFORCE_EDDINGTON_LIMIT # put a hard limit on the maximum accretion rate
#BH_THERMALFEEDBACK      # couple a fraction of the BH luminosity into surrounding 
#BH_DRAG                 # Drag on black-holes due to accretion
#SWALLOWGAS              # Enables stochastic accretion of gas particles consistent with growth rate of hole
#DRAINGAS=1              # non-stochastic/smooth draining of accreted mass (1: smooth draining, 2: smooth draing + rho for Bondi rate based on primary cell)
#REPOSITION_ON_POTMIN    # repositions hole on potential minimum (requires EVALPOTENTIAL)
#BH_COUNTPROGS	         # carries a counter for each BH that gives the total number of seeds that merged into it
#BH_DYN_FRICTION
#BH_USE_GASVEL_IN_BONDI  # only when this is enabled, the surrounding gas velocity is used in addition to the sounds speed in the Bondi rate 
#CSND_FRAC_BH_MERGE=0.5  # Relative velocity fraction (in units of soundspeed) for merging blackholes, default=0.5
#BH_CAPTURE_RATE         # Alternative to the Bondi prescription
#BONDI_DISK_VORTICITY         # Needs to be activated together with BONDI or with BH_CAPTURE_RATE
#MASSIVE_SEEDS           # BH seeds assigned large dynamical mass, such that ideally no repositioning is needed anymore
#-------------------------------------- AGN-Bubble feedback
#BH_BUBBLES              # calculate bubble energy directly from the black hole accretion rate
#UNIFIED_FEEDBACK        # activates BH_THERMALFEEDBACK at high Mdot and BH_BUBBLES FEEDBACK al low Mdot

#---------------------------------------- Passive Tracers
#TRACER_FIELD
#TRACER_PARTICLE=2                 # advect massless tracer particles of type TRACER_PARTICLE with velocity field
#TRACER_MC=1 #=2                   # Monte Carlo tracer particles (=1 to enable, >=2 to output as that partType)
#GENERATE_TRACER_PARTICLE_IN_ICS   # add tracer particles at positions of cell vertices in ICs
#TRACER_MC_DEBUG
#TRACER_MC_NOEXCHANGE

#TRACER_PART_NUM_FLUID_QUANTITIES=6     # number of fluid quantities to be stored for velocity tracers - must match the value given to TRACER_PART_STORE_WHAT
#TRACER_PART_STORE_WHAT=1+2+4+8+16+32   # bit mask for quantities to store (see allvars.h for bitmask)

#TRACER_MC_NUM_FLUID_QUANTITIES=9       # number of fluid quantities to be stored for MC tracers - must match the value gien to TRACER_MC_STORE_WHAT
#TRACER_MC_STORE_WHAT=1+2+4+8+16+32+64+128+256     # bit mask for quantities to store (see allvars.h for bitmask)

#TRACER_MC_SKIPLOAD=3                   # skip reading this particle type when reading initial conditions from a snapshot file
#FOF_DISABLE_SNAP_REWRITE               # do not rewrite a snapshot file when RestartFlag==3


#-------------------------------------------- Things for special behaviour
#TILE_ICS
#USE_RANDOM_GENERATOR
#MULTIPLE_RESTARTS
#TOLERATE_WRITE_ERROR
#OPTIMIZE_MEMORY_USAGE                   #optimize for memory, not for speed. Note: this is dangerous for high dynamic range simulations with mixed precision, since some position variables are singles instead of doubles
#SUBBOX_SNAPSHOTS
#INDIVIDUAL_GRAVITY_SOFTENING
#PROCESS_TIMES_OF_OUTPUTLIST
#EXTENDED_GHOST_SEARCH           # This extends the ghost search to the full 3x3 domain instead of the principal domain
#ALTERNATIVE_GHOST_SEARCH        # This switches on the "old" routines that find the ghost neighbours

VORONOI_DYNAMIC_UPDATE          # keeps track of mesh connectivity, which speeds up mesh construction
#COFFEE_PROBLEM
#NOH_PROBLEM
#SHIFT_BY_HALF_BOX
#DISABLE_VELOCITY_CSND_SLOPE_LIMITING
NO_FANCY_MPI_CONSTRUCTS
NO_MPI_IN_PLACE
NO_ISEND_IRECV_IN_DOMAIN
FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
#MPISENDRECV_SIZELIMIT=100
#MPI_HYPERCUBE_ALLGATHERV     # some MPI-libraries may use quite a bit of internal storage for MPI_Allgatherv. This uses hypercubes instead as a work-around
#MPISENDRECV_CHECKSUM
#NOTREERND
#ENLARGE_DYNAMIC_RANGE_IN_TIME   # This extends the dynamic range of the integer timeline from 32 to 64 bit
#INCLUDE_VMAX_IN_TREEBASED_TISTEP_OPENING 
#NOSTOP_WHEN_BELOW_MINTIMESTEP
#NOTYPEPREFIX_FFTW
#DO_NOT_CREATE_STAR_PARTICLES
#ALLOWEXTRAPARAMS
#NEW_RATES                     # switches in updated cooling rates from Naoki
#RADIATIVE_RATES               # used in non-equilibrium chemistry model
#FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES  # this can be used to load SPH ICs that contain identical particle coordinates
#VEL_POWERSPEC                 # compiles in a code module that allows via restart-flag 7 the calculation of a gas velocity power spectrum of a snapshot
#ADJ_BOX_POWERSPEC         # compiles in a code module that allows via restart-flag 7 the calculation of gas power spectra of a snapshot with an adjustable box (user defined center and size)
#ALTERNATIVE_FFT
#DISABLE_OPTIMIZE_DOMAIN_MAPPING

#--------------------------------------- Output/Input options
#REDUCE_FLUSH
#OUTPUT_EVERY_STEP
#GODUNOV_STATS
#OUTPUT_PRESSURE_GRADIENT
#OUTPUT_DENSITY_GRADIENT
#OUTPUT_VELOCITY_GRADIENT
#OUTPUT_VERTEX_VELOCITY
#OUTPUT_VERTEX_VELOCITY_DIVERGENCE  # requires CALCULATE_VERTEX_VELOCITY_DIVERGENCE
#OUTPUT_CENTER_OF_MASS
#OUTPUT_SURFACE_AREA
#OUTPUT_PRESSURE
#OUTPUTPOTENTIAL
#OUTPUTACCELERATION
#OUTPUTTIMESTEP
HAVE_HDF5                     # needed when HDF5 I/O support is desired
#PARAMS_IN_SNAP                # add the compiler flags and parameter file values to every snapshot file (requires HAVE_HDF5)
#HDF5_FILTERS                  # activate snapshot compression and checksum for HDF5 output
#OUTPUTCOOLRATE                # outputs cooling rate, and conduction rate if enabled
#OUTPUT_DIVVEL                 # output  velocity divergence
#OUTPUT_CURLVEL                 # output  velocity curl
#OUTPUT_COOLHEAT               # output actual energy loss/gain in cooling/heating routine
#OUTPUT_VORTICITY
#OUTPUT_CELL_SPIN              
#MEASURE_DISSIPATION_RATE      # measures and outputs dissipation rate. Note: requires USE_ENTROPY_FOR_COLD_FLOWS, even though it will then always use the thermal energy update
#OUTPUT_MACHNUM                # output maximum mach number of a cell

#--------------------------------------- Testing and Debugging options
#DEBUG                         # enables core-dumps 
#DEBUG_ENABLE_FPU_EXCEPTIONS   # tries to enable FPU exceptions
#VERBOSE                       # reports readjustments of buffer sizes
#HOST_MEMORY_REPORTING         # reports after start-up the available system memory by analyzing /proc/meminfo
#FORCETEST                     # calculates for given fraction of particles direct summation forces to check accuracy of tree force

#--------------------------------------- Static Disk Potential
#DISK_POTENTIAL
#DISK_MASS_M0=1.0
#DISK_SCALE_R0=1.0

#--------------------------------------- Static NFW Potential
#STATICNFW
#NFW_C=12
#NFW_M200=100.0
#NFW_Eps=0.01
#NFW_DARKFRACTION=0.87

#--------------------------------------- Static Isothermal Sphere Potential
#STATICISO                                                                                                                                       
#ISO_M200=100.0                                                                                                                                  
#ISO_R200=160.0                                                                                                                                  
#ISO_Eps=0.1                                                                                                                                     
#ISO_FRACTION=0.9  

#--------------------------------------- Static Hernquist Potential
#STATICHQ
#HQ_M200=186.015773
#HQ_C=10.0
#HQ_DARKFRACTION=0.9

#--------------------------------------- Growing Disk Potential
#GROWING_DISK_POTENTIAL

#--------------------------------------- Dark energy
#DARKENERGY # Enables Dark Energy
#TIMEDEPDE  # read w(z) from a DE file
#RESCALEVINI # rescale v_ini in read_ic / read_ic_cluster
#EXTERNALHUBBLE # reads the hubble function from the DE file
#TIMEDEPGRAV # resacles H and G according to DE model
#DARKENERGY_DEBUG # enable writing of drift/kick table

#--------------------------------------- Glass making/ 2nd-order initial conditions / Initial conditions options
#SECOND_ORDER_ICS
#LONGIDS
#INHOMOG_GASDISTR_HINT         # if the gas is distributed very different from collisionless particles, this can helps to avoid problems in the domain decomposition
#GENERATE_GAS_IN_ICS
#SPLIT_PARTICLE_TYPE=4+8

#-------------------------------------- Simple turbulence test
#VS_TURB
#POWERSPEC_GRID=128
#GAMMA=1.01

#AB_TURB

#--------------------------------------- Planets/Materials (Robert)
#MATERIALS                   # master switch for all extension for planetary physics 
#NUMBER_OF_MATERIALS=1
#READ_LEGACY_ICS

#--------------------------------------- Degenerate Equation of State
#EOS_DEGENERATE
#EOS_COULOMB_CORRECTIONS
#EOS_NSPECIES=3
#RELAXOBJECT
#PASSIVE_SCALARS=3

#-------------------------------------- Nuclear Network
#NUCLEAR_NETWORK
#REACLIB1

#-------------------------------------- Radiative transfer options
#RT_COOLING_PHOTOHEATING
#RT_ADVECT                   # enable advection of radiation field
#RT_CGMETHOD                   # enables CG method solution of the RT advection
#RT_SLOWLIGHT                # enable slow light approximation
#RT_N_DIR=2                 # track this number of locally brightest sources (one of them is diffuse field)
#RT_COMBINE_N_DIR_IN_OUTPUT  # writes only a single summed photon/photon-density field into output files
#RT_ALLOW_ABSORBING_CELLS    # if this is set, all cells with ID>=1000000000 will absorb radiation
#RT_SPREAD_SOURCE
#RT_STELLAR_SOURCES
#RT_HEALPIX_NSIDE=1          # if this is set, a discretization of the solid angle is used instead of brightest source selection
#RT_INCLUDE_HE
#SOURCE_PERIODIC

#DO_NOT_MOVE_GAS
#HYDROGEN_ONLY

#-------------------------------------- TG options
#TGSET
#TG_READ_UTHERM
#PRIMCHEM
#PRIMCHEM_OUTRATES
#TGCOOL
#RELATIVE_VELOCITY

#-------------------------------------- Calculate Hessian Matrix
#SECOND_DERIVATIVES
#SLOPE_LIMIT_HESSIANS
#RECONSTRUCT_GRADIENTS

#-------------------------------------- Navier-Stokes Terms
#GLOBAL_VISCOSITY 	   #needs dynamic and bulk coefficients
#USE_KINEMATIC_VISCOSITY    #needs only one input parameter
#ALPHA_VISCOSITY            #for accretion disks
#LOCAL_VISCOSITY=1          #=1 Sutherland viscosity/ =2 Spitzer viscosity
#THERMAL_CONDUCTION
#TRACER_DIFFUSION           #requires TRACER_FIELD switched on

#-------------------------------------- Circumstellar Disks
#CIRCUMSTELLAR              #Master switch
#CENTRAL_MASS_POTENTIAL     #Point-mass potential
#STAR_PLANET_POTENTIAL     #Fixed star-planet circular orbit
#LOCALLY_ISOTHERM_DISK      #Isothermal Equation of state at each radii.

#-------------------------------------- Special Boundaries within domain
#SPECIAL_BOUNDARY          #Main Switch
#COAXIAL_BOUNDARIES              #e.g. Couette flow-type boundaries

#-------------------------------------- Windtunnel
#WINDTUNNEL
#BOUNDARY_INFLOW_OUTFLOW_MIN=30000000  # defines the ID range describing inflow/outflow nozzle of wind-tunnel
#BOUNDARY_INFLOW_OUTFLOW_MAX=40000000
#WINDTUNNEL_COORD=0                    # sets the coordinate in which the wind blows (0,1,2 for x,y,z)
#WINDTUNNEL_EXTERNAL_SOURCE

#-------------------------------------- GFM - Galaxy Formation Module
#GFM                                    #master switch
#GFM_STELLAR_EVOLUTION=0                #stellar evolution: 0->default, 1->no mass loss (note: beta value changes + MassMetallicity & MassMetals inconsistent internally with cell dynamical mass) 2->call only test routine 
#GFM_STELLAR_FEEDBACK                   #SNIa and AGB feedback
#GFM_PREENRICH                          #pre enrich gas at given redshift
#GFM_WINDS=0                            #decoupled winds: =0 stochastically, >0 spawn wind particle from star-forming gas every GFM_WINDS-th step wind is considered
#GFM_WINDS_VARIABLE=0                   #0->scale winds with halo mass, requires FoF, 1->sigma winds
#GFM_BIPOLAR_WINDS=1                    #bipolar winds: 0->default, 1->relative to motion of FOF group
#GFM_COOLING_METAL=0                    #metal line cooling: 0->scale with metallicity, 1->interpolate metallicity
#GFM_OUTPUT_MASK=1+2+4+8+16             #which fields to output ()

#-------------------------------------- FM - Star formation and feedback module
#FM_SFR                                #turns on star formation (needs USE_SFR)
#FM_STELLAR_FEEDBACK                   #turns on stellar feedback
#DELAYED_COOLING                       #turns on delayed cooling model (Stinson et al. 2006)
#EXPLICIT_COOLING                      #switch to a 2-nd order explicit method for cooling if (u^{n+1} - u_{n}) < tol * u^{n}
#COMPUTE_SFR_FROM_H2                   #links the SFR to the H2 gas fraction
#TEST_COOLING_METAL                    #call only cooling test routine (save cooling function with metal cooling for solar metallicity)
#OUTPUT_STELLAR_FEEDBACK               #outputs SNII number, feedback energy and mass released for stellar particles and log files for feedback (requires GFM_STELLAR_EVOLUTION)
#OUTPUT_MOLECULAR_FRACTION             #outputs the H2 gas fraction (requires COMPUTE_SFR_FROM_H2 switched on)
#OUTPUT_OPTICAL_DEPTH                  #outputs the gas optical depth (requires COMPUTE_SFR_FROM_H2 switched on)

