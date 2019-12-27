/*
 * allvars.h derived
 * dnelson
 */

#ifndef AREPO_RT_ALLVARS_H
#define AREPO_RT_ALLVARS_H

/* variable types (including dtypes.h) */

#define  MAX_FLOAT_NUMBER 1e37
#define  MIN_FLOAT_NUMBER 1e-37
#define  MAX_DOUBLE_NUMBER 1e306
#define  MIN_DOUBLE_NUMBER 1e-306

#define  NTYPES 6
#define  NSOFTTYPES NTYPES
#define  NSOFTTYPES_HYDRO 0

#ifdef DOUBLEPRECISION
#if (DOUBLEPRECISION==2)
#define  MAX_REAL_NUMBER  MAX_FLOAT_NUMBER
#define  MIN_REAL_NUMBER  MIN_FLOAT_NUMBER
#else
#define  MAX_REAL_NUMBER  MAX_DOUBLE_NUMBER
#define  MIN_REAL_NUMBER  MIN_DOUBLE_NUMBER
#endif
#else
#define  MAX_REAL_NUMBER  MAX_FLOAT_NUMBER
#define  MIN_REAL_NUMBER  MIN_FLOAT_NUMBER
#endif

/* memory management */

#ifndef DISABLE_MEMORY_MANAGER
#define  mymalloc(x, y)            mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  mymalloc_movable(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__)

#define  myrealloc(x, y)           myrealloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  myrealloc_movable(x, y)   myrealloc_movable_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)

#define  myfree(x)                 myfree_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
#define  myfree_movable(x)         myfree_movable_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)

#define  report_memory_usage(x, y) report_detailed_memory_usage_of_largest_task(x, y, __FUNCTION__, __FILE__, __LINE__)

#else
#define  mymalloc(x, y)            malloc(y)
#define  mymalloc_movable(x, y, z) malloc(z)

#define  myrealloc(x, y)           realloc(x, y)
#define  myrealloc_movable(x, y)   realloc(x, y)

#define  myfree(x)                 free(x)
#define  myfree_movable(x)         free(x)

#define  report_memory_usage(x, y) printf("Memory manager disabled.\n")
#endif

/* magic constants */

#define MAXLEN_OUTPUTLIST    1100
#define TIMEBINS             29
#define GRAVCOSTLEVELS       6
#define CPU_LAST             43

#define NUMBER_OF_MEASUREMENTS_TO_RECORD  6

#define FACT1                0.366025403785    /* FACT1 = 0.5 * (sqrt(3)-1) */
#define MAXLEN_PATH          256

typedef int integertime;

/* sph kernel interpolation */
#define  NUMDIMS 3              /**< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470 /**< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786 /**< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */

/* box and periodic */

extern MyDouble boxSize, boxHalf; // NOTE: this is maybe not the boxSize global used by Arepo functions since this is
                                  // "C++ mangled" (whether or not a function is inside an extern "C" determines which it sees
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf

#ifndef REFLECTIVE_X
#define NGB_PERIODIC_LONG_X(x) (xtmp = fabs(x), (xtmp > boxHalf_X) ? (boxSize_X - xtmp) : xtmp)
#define NEAREST_X(x) (xtmp = (x), (xtmp > boxHalf_X) ? (xtmp - boxSize_X) : ((xtmp < -boxHalf_X) ? (xtmp + boxSize_X) : (xtmp)))
#define WRAP_X(x) (xtmp = (x), (xtmp > boxSize_X) ? (xtmp - boxSize_X) : ((xtmp < 0) ? (xtmp + boxSize_X) : (xtmp)))
#else /* #ifndef REFLECTIVE_X */
#define NGB_PERIODIC_LONG_X(x) fabs(x)
#define NEAREST_X(x) (x)
#define WRAP_X(x) (x)
#endif /* #ifndef REFLECTIVE_X #else */

#ifndef REFLECTIVE_Y
#define NGB_PERIODIC_LONG_Y(x) (ytmp = fabs(x), (ytmp > boxHalf_Y) ? (boxSize_Y - ytmp) : ytmp)
#define NEAREST_Y(x) (ytmp = (x), (ytmp > boxHalf_Y) ? (ytmp - boxSize_Y) : ((ytmp < -boxHalf_Y) ? (ytmp + boxSize_Y) : (ytmp)))
#define WRAP_Y(x) (ytmp = (x), (ytmp > boxSize_Y) ? (ytmp - boxSize_Y) : ((ytmp < 0) ? (ytmp + boxSize_Y) : (ytmp)))
#else /* #ifndef REFLECTIVE_Y */
#define NGB_PERIODIC_LONG_Y(x) fabs(x)
#define NEAREST_Y(x) (x)
#define WRAP_Y(x) (x)
#endif /* #ifndef REFLECTIVE_Y #else */

#ifndef REFLECTIVE_Z
#define NGB_PERIODIC_LONG_Z(x) (ztmp = fabs(x), (ztmp > boxHalf_Z) ? (boxSize_Z - ztmp) : ztmp)
#define NEAREST_Z(x) (ztmp = (x), (ztmp > boxHalf_Z) ? (ztmp - boxSize_Z) : ((ztmp < -boxHalf_Z) ? (ztmp + boxSize_Z) : (ztmp)))
#define WRAP_Z(x) (ztmp = (x), (ztmp > boxSize_Z) ? (ztmp - boxSize_Z) : ((ztmp < 0) ? (ztmp + boxSize_Z) : (ztmp)))
#else /* #ifndef REFLECTIVE_Z */
#define NGB_PERIODIC_LONG_Z(x) fabs(x)
#define NEAREST_Z(x) (x)
#define WRAP_Z(x) (x)
#endif /* #ifndef REFLECTIVE_Z #else */

#define ALLOC_TOLERANCE 0.1

/* file pointers */

extern FILE *FdInfo,            /**< file handle for info.txt log-file. */
 *FdEnergy,                     /**< file handle for energy.txt log-file. */
 *FdTimings,                    /**< file handle for timings.txt log-file. */
 *FdBalance,                    /**< file handle for balance.txt log-file. */
 *FdTimebin,                    /**< file handle for timebins.txt log-file. */
 *FdDomain,                     /**< file handle for domain.txt log-file. */
 *FdMemory,                     /**< file handle for memory.txt log-file. */
 *FdCPU;                        /**< file handle for cpu.txt log-file. */


/*********************************************************/
/*  Global variables                                     */
/*********************************************************/

extern int ThisTask, NTask, PTask;
extern int NumPart, NumGas;
extern int RestartFlag, RestartSnapNum, WriteMiscFiles;
extern char ParameterFile[MAXLEN_PATH];
extern int TimeBinActive[TIMEBINS];

extern struct global_data_all_processes
{
  // ArepoVTK note: this has to match with the compiled Arepo version to read the parameterfile ok
	
  long long TotNumPart;		/**<  total particle numbers (global value) */
  long long TotNumGas;		/**<  total gas particle number (global value) */

  int MaxPart;
  int MaxPartSph;
  
#ifdef COOLING
  char TreecoolFile[MAXLEN_PATH];
#endif

  double TotGravCost;
  double MeanVolume;
  int    MultipleDomains;
  double TopNodeFactor;
  int ICFormat;
  int SnapFormat;
  int NumFilesPerSnapshot;
  int NumFilesWrittenInParallel;
  double TreeAllocFactor;
  double TopNodeAllocFactor;
  double NgbTreeAllocFactor;
  int MaxMemSize;

  /* some SPH parameters */

  int DesNumNgb;		/**< Desired number of SPH neighbours */
  double TotCountReducedFluxes;
  double TotCountFluxes;
  double DtDisplacement;
  double MaxNumNgbDeviation;	/**< Maximum allowed deviation neighbour number */
  double InitGasTemp;		/**< may be used to set the temperature in the IC's */
  double InitGasU;		/**< the same, but converted to thermal energy per unit mass */
  double MinGasTemp;		/**< may be used to set a floor for the gas temperature */
  double MinEgySpec;		/**< the minimum allowed temperature expressed as energy per unit mass */

  double MinimumDensityOnStartUp;
  double GasSoftFactor;
  double LimitUBelowThisDensity;
  double LimitUBelowCertainDensityToThisValue;

  /* some force counters  */

  long long TotNumOfForces;	/**< counts total number of force computations  */

  double cf_atime, cf_a2inv, cf_a3inv, cf_afac1, cf_afac2, cf_afac3, cf_hubble_a, cf_time_hubble_a, cf_redshift;
	double cf_H;
  double cf_Hrate;

  /* system of units  */
  double UnitTime_in_s,		/**< factor to convert internal time unit to seconds/h */
    UnitMass_in_g,		/**< factor to convert internal mass unit to grams/h */
    UnitVelocity_in_cm_per_s,	/**< factor to convert internal velocity unit to cm/sec */
    UnitLength_in_cm,		/**< factor to convert internal length unit to cm/h */
    UnitPressure_in_cgs,	/**< factor to convert internal pressure unit to cgs units (little 'h' still
				   around!) */
    UnitDensity_in_cgs,		/**< factor to convert internal length unit to g/cm^3*h^2 */
    UnitCoolingRate_in_cgs,	/**< factor to convert internal cooling rate to cgs units */
    UnitEnergy_in_cgs,		/**< factor to convert internal energy to cgs units */
    UnitTime_in_Megayears,	/**< factor to convert internal time to megayears/h */
    GravityConstantInternal,	/**< If set to zero in the parameterfile, the internal value of the
				   gravitational constant is set to the Newtonian value based on the system of
				   units specified. Otherwise the value provided is taken as internal gravity
				   constant G. */
    G;				/**< Gravity-constant in internal units */

  /* Cosmology */

  double Hubble;		/**< Hubble-constant in internal units */
  double Omega0,		/**< matter density in units of the critical density (at z=0) */
    OmegaLambda,		/**< vaccum energy density relative to crictical density (at z=0) */
    OmegaBaryon,		/**< baryon density in units of the critical density (at z=0) */
    HubbleParam;		/**< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute
				 * physical values for cooling physics
				 */

  double BoxSize;		/**< Boxsize in case periodic boundary conditions are used */

  int ComovingIntegrationOn;	/**< flags that comoving integration is enabled */
  int PeriodicBoundariesOn;	/**< flags that periodic boundaries are enabled */
  int ResubmitOn;		/**< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;	/**< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative
				   criterion */
  int TypeOfTimestepCriterion;	/**< gives type of timestep criterion (only 0 supported right now - unlike
				   gadget-1.1) */
  int OutputListOn;		/**< flags that output times are listed in a specified file */
  int CoolingOn;		/**< flags that cooling is enabled */
  int StarformationOn;		/**< flags that star formation is enabled */

  int NParameters;

  int LowestActiveTimeBin;
  int HighestActiveTimeBin;
  int LowestOccupiedTimeBin;
  int HighestOccupiedTimeBin;
  int LowestOccupiedGravTimeBin;
  int HighestOccupiedGravTimeBin;
  int HighestSynchronizedTimeBin;
  int SmallestTimeBinWithDomainDecomposition;
  double ActivePartFracForNewDomainDecomp;

  int SnapshotFileCount;	/**< number of snapshot that is written next */

  double TimeBetSnapshot,	/**< simulation time interval between snapshot files */
    TimeOfFirstSnapshot,	/**< simulation time of first snapshot files */
    CpuTimeBetRestartFile,	/**< cpu-time between regularly generated restart files */
    TimeLastRestartFile,	/**< cpu-time when last restart-file was written */
    TimeBetStatistics,		/**< simulation time interval between computations of energy statistics */
    TimeLastStatistics;		/**< simulation time when the energy statistics was computed the last time */
  int NumCurrentTiStep;		/**< counts the number of system steps taken up to this point */

  /* Current time of the simulation, global step, and end of simulation */

  double Time,			/**< current time of the simulation */
    TimeBegin,			/**< time of initial conditions of the simulation */
    TimeStep,			/**< difference between current times of previous and current timestep */
    TimeMax;			/**< marks the point of time until the simulation is to be evolved */

  /* variables for organizing discrete timeline */

  double Timebase_interval;	/**< factor to convert from floating point time interval to integer timeline */
  integertime Ti_Current;	/**< current time on integer timeline */
  integertime Previous_Ti_Current;
  integertime Ti_nextoutput;	/**< next output time on integer timeline */
  integertime Ti_lastoutput;

  integertime Ti_begstep[TIMEBINS];    /**< marks start of current step of each timebin on integer timeline */

  long long GlobalNSynchronizedHydro;
  long long GlobalNSynchronizedGravity;

  int LevelToTimeBin[GRAVCOSTLEVELS];
  int LevelHasBeenMeasured[GRAVCOSTLEVELS];

  /* variables that keep track of cumulative CPU consumption */
  double TimeLimitCPU;
  double CPU_Sum[CPU_LAST];

  /* tree code opening criterion */
  double ErrTolTheta;
  double ErrTolForceAcc;

  /* adjusts accuracy of time-integration */
  double ErrTolIntAccuracy;
  double MinSizeTimestep,
    MaxSizeTimestep;

  double IsoSoundSpeed;
  double CourantFac;

#ifdef REGULARIZE_MESH_FACE_ANGLE
  double CellMaxAngleFactor;
#else
  double CellShapingFactor;
#endif
  double CellShapingSpeed;

  int CPU_TimeBinCountMeasurements[TIMEBINS];
  double CPU_TimeBinMeasurements[TIMEBINS][NUMBER_OF_MEASUREMENTS_TO_RECORD];

  int SofteningTypeOfPartType[NTYPES];

  double SofteningComoving[NSOFTTYPES];
  double SofteningMaxPhys[NSOFTTYPES];

  double
      SofteningTable[NSOFTTYPES + NSOFTTYPES_HYDRO];
  double ForceSoftening[NSOFTTYPES + NSOFTTYPES_HYDRO + 1];
  double MassTable[NTYPES];

  /* some filenames */
  char InitCondFile[MAXLEN_PATH],
    OutputDir[MAXLEN_PATH],
    SnapshotFileBase[MAXLEN_PATH],
    ResubmitCommand[MAXLEN_PATH],
    OutputListFilename[MAXLEN_PATH];

  /** table with desired output times */
  double OutputListTimes[MAXLEN_OUTPUTLIST];
  char OutputListFlag[MAXLEN_OUTPUTLIST];
  int OutputListLength;		/**< number of times stored in table of desired output times */

  MyIDType MaxID;

#ifdef SPECIAL_BOUNDARY
  double BoundaryLayerScaleFactor;
  double SpecialBoundarySpeed;
  int SpecialBoundaryMotion;
  int SpecialBoundaryType;
  double OutflowPressure;
#endif

  double GlobalDisplacementVector[3];
}
All;

extern struct particle_data
{
  MyDouble Pos[3]; // __attribute__((__aligned__(16)));
  MyDouble Mass;
  MyFloat  Vel[3]; // __attribute__((__aligned__(16)));
  MyFloat  GravAccel[3];
  MyIDType ID;

  integertime Ti_Current; /*!< current time on integer timeline */

  float OldAcc; // ArepoVTK: used to store ElectronAbundance (Ne)
	
  float GravCost[GRAVCOSTLEVELS];	/**< weight factors used for balancing the work-load */
  unsigned char Type;		/**< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
  unsigned char SofteningType;
  signed char TimeBinGrav;
  signed char TimeBinHydro;
}
 *P,				/**< holds particle data on local processor */
 *DomainPartBuf;		/**< buffer for particle data used in domain decomposition */

extern struct sph_particle_data
{
  /* conserved variables */
  MyFloat Energy;
  MyFloat Momentum[3];
  MyFloat Volume;
  MyFloat OldMass; // ArepoVTK: used to store Entropy (code units)

  /* primitive variables */
  MyFloat Density;
  MyFloat Pressure;		/**< current pressure */
  MySingle Utherm;

  /* variables for mesh  */
  MyDouble Center[3];		/* center of mass of cell */
  MySingle VelVertex[3];		/* current vertex velocity (primitive variable) */ /* ArepoVTK: used to store {Bmag,ShockDEDT,empty} */
  MySingle MaxDelaunayRadius;
  MySingle Hsml;		        /* auxiliary search radius for points around a delaunay triangle */
  MySingle SurfaceArea;
  MySingle ActiveArea;

#ifdef METALS
  MyFloat MassMetallicity;
#endif

#ifdef COOLING
  MyFloat Ne;
#endif

#ifdef USE_SFR
  MySingle Sfr;
#endif

  struct grad_data Grad;

#ifdef SECOND_DERIVATIVES
  struct hessian_data Hessian;
#endif

  int first_connection;
  int last_connection;

#ifdef SPECIAL_BOUNDARY
  MyFloat MinDistBoundaryCell;
#endif
  double TimeLastPrimUpdate;

}
 *SphP,				/**< holds SPH particle data on local processor */
 *DomainSphBuf;			/**< buffer for SPH particle data in domain decomposition */

extern int NTopnodes, NTopleaves;
 
extern struct NODE
{
  union
  {
    int suns[8];
    struct
    {
      MyDouble s[3];
      MyDouble mass;
      int sibling;		
      int nextnode;
      int father;
      unsigned char maxsofttype; // #if(NSOFTTYPES > 1)
    }
    d;
  }
  u;

  MyDouble center[3];           /**< geometrical center of node */
  MyFloat len;                 /**< sidelength of treenode */
}
 *Nodes;	
 
extern int *Nextnode;		
extern int *Father;		 
 
/** Variables for neighbor tree */
extern int Ngb_MaxPart;
extern int Ngb_NumNodes;
extern int Ngb_MaxNodes;
extern int Ngb_FirstNonTopLevelNode;
extern int Ngb_NextFreeNode;
extern int *Ngb_Father;
extern int *Ngb_Marker;
extern int Ngb_MarkerValue;

extern int *Ngb_DomainNodeIndex;
extern int *DomainListOfLocalTopleaves;
extern int *DomainNLocalTopleave;
extern int *DomainFirstLocTopleave;
extern int *Ngb_Nextnode;

/** The ngb-tree data structure
 */
extern struct NgbNODE
{
  union
  {
    int suns[8];                /**< temporary pointers to daughter nodes */
    struct
    {
      int sibling;
      int nextnode;
      MyNgbTreeFloat range_min[3];
      MyNgbTreeFloat range_max[3];
    } d;
  } u;

  MyNgbTreeFloat vertex_vmin[3];
  MyNgbTreeFloat vertex_vmax[3];

  int father;

  integertime Ti_Current;
}
 *Ngb_Nodes; 
 
/*
 * proto.h derived
 */

extern "C" { 
void begrun1(void);
void endrun(int);
int init(void); 
void allocate_memory(void);
void set_cosmo_factors_for_current_time(void);
void read_ic(const char *fname, int); 
void determine_compute_nodes(void);
 
void domain_Decomposition(void);
void set_softenings(void);
int ngb_treebuild(int npart);
void ngb_treeallocate(void);
void setup_smoothinglengths(void);
void density(void);
void reconstruct_timebins(void);
 
void create_mesh(void); 
void voronoi_init_connectivity(tessellation * T);

void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int linenr);
void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file,int line);

void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line);
void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line);

void myfree_fullinfo(void *p, const char *func, const char *file, int line);
void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line);

void mymalloc_init(void);
void dump_memory_table(void);

void open_logfiles(void);
void close_logfiles(void);
void write_voronoi_mesh(tessellation * T, char *fname, int writeTask, int lastTask);

int image_get_next_tetra(tessellation * T, int tt, point * ppstart, point * ppend, int *nexttetra,
                         point * ppexit, int *previous_tetra);

void calc_mpz_determinant(mpz_t det, mpz_t ax, mpz_t ay, mpz_t az, mpz_t bx, mpz_t by, mpz_t bz, mpz_t cx, mpz_t cy, mpz_t cz);

}
 
#endif //AREPO_RT_ALLVARS_H
