/*
 * allvars.h derived (July 2012)
 * dnelson
 */

#ifndef AREPO_RT_ALLVARS_H
#define AREPO_RT_ALLVARS_H

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

#define MAX_REAL_NUMBER      1e37
#define MAXLEN_OUTPUTLIST    1100
#define TIMEBINS             60
#define GRAVCOSTLEVELS       6
#define CPU_PARTS            51

#define NUMBER_OF_MEASUREMENTS_TO_RECORD  6

#define FACT1                0.366025403785    /* FACT1 = 0.5 * (sqrt(3)-1) */
#define SUNRISE              4711
#define MAXLEN_PATH          256

typedef int integertime;

/* box and periodic */

extern MyDouble boxSize, boxHalf;
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf

#ifdef PERIODIC
#define NGB_PERIODIC_LONG_X(x) (xtmp=fabs(x),(xtmp>boxHalf_X)?(boxSize_X-xtmp):xtmp)
#define NGB_PERIODIC_LONG_Y(x) (ytmp=fabs(x),(ytmp>boxHalf_Y)?(boxSize_Y-ytmp):ytmp)
#define NGB_PERIODIC_LONG_Z(x) (ztmp=fabs(x),(ztmp>boxHalf_Z)?(boxSize_Z-ztmp):ztmp)

#define NEAREST_X(x) (xtmp=(x),(xtmp>boxHalf_X)?(xtmp-boxSize_X):( (xtmp< -boxHalf_X)?(xtmp+boxSize_X):(xtmp) ) )
#define NEAREST_Y(x) (ytmp=(x),(ytmp>boxHalf_Y)?(ytmp-boxSize_Y):( (ytmp< -boxHalf_Y)?(ytmp+boxSize_Y):(ytmp) ) )
#define NEAREST_Z(x) (ztmp=(x),(ztmp>boxHalf_Z)?(ztmp-boxSize_Z):( (ztmp< -boxHalf_Z)?(ztmp+boxSize_Z):(ztmp) ) )
#else

#define NGB_PERIODIC_LONG_X(x) fabs(x)
#define NGB_PERIODIC_LONG_Y(x) fabs(x)
#define NGB_PERIODIC_LONG_Z(x) fabs(x)

#define NEAREST_X(x) (x)
#define NEAREST_Y(x) (x)
#define NEAREST_Z(x) (x)
#endif

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

extern struct global_data_all_processes
{
  // ArepoVTK note: this has to match with the compiled Arepo version to read the parameterfile ok
	
  long long TotNumPart;		/**<  total particle numbers (global value) */
  long long TotNumGas;		/**<  total gas particle number (global value) */

  int MaxPart;
  int MaxPartSph;
  
#ifdef COOLING
  char TreecoolFile[255];
#endif

#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS) 
  double TargetGasMass;
  int RefinementCriterion;
  int DerefinementCriterion;
#endif

  double TotGravCost;
  double MeanVolume;
  int    MultipleDomains;
  double TopNodeFactor;
  int ICFormat;
  int SnapFormat;
  int NumFilesPerSnapshot;
  int NumFilesWrittenInParallel;
  int BufferSize;
  int BunchSize;
  double TreeAllocFactor;
  double TopNodeAllocFactor;
  double NgbTreeAllocFactor;
  int MaxMemSize;

  /* some SPH parameters */

  int DesNumNgb;		/**< Desired number of SPH neighbours */
  double TotCountReducedFluxes;
  double TotCountFluxes;
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

#ifdef SUBBOX_SNAPSHOTS
  double SubboxXmin, SubboxXmax, SubboxYmin, SubboxYmax, SubboxZmin, SubboxZmax;
  double SubboxMinTime, SubboxMaxTime; 
  int SubboxSyncCounter;
  int SubboxSyncModulo;
  int SubboxNumFilesPerSnapshot;      
  int SubboxNumFilesWrittenInParallel;
#endif

  double cf_atime, cf_a2inv, cf_a3inv, cf_afac1, cf_afac2, cf_afac3, cf_hubble_a;

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

  int LowestActiveTimeBin;
  int HighestActiveTimeBin;
  int HighestOccupiedTimeBin;
  int SmallestTimeBinWithDomainDecomposition;
  int MaxTimeBinsWithoutDomainDecomposition;

  int SnapshotFileCount;	/**< number of snapshot that is written next */
#ifdef SUBBOX_SNAPSHOTS
  int SubboxSnapshotFileCount;  /**< number of subbox snapshot that is written next */
#endif
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

  int LevelToTimeBin[GRAVCOSTLEVELS];
  int LevelHasBeenMeasured[GRAVCOSTLEVELS];

  /* variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;
  double CPU_Sum[CPU_PARTS];	/**< sums wallclock time/CPU consumption in whole run */

  /* tree code opening criterion */

  double ErrTolTheta;		/**< BH tree opening angle */
  double ErrTolForceAcc;	/**< parameter for relative opening criterion in tree walk */

  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/**< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
				   timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep,	/**< minimum allowed timestep. Normally, the simulation terminates if the
				   timestep determined by the timestep criteria falls below this limit. */
    MaxSizeTimestep;		/**< maximum allowed timestep */

  double MaxRMSDisplacementFac;	
  double IsoSoundSpeed;
  double CourantFac;		/**< SPH-Courant factor */

#ifdef REGULARIZE_MESH_FACE_ANGLE
  double CellMaxAngleFactor;
#else
  double CellShapingFactor;
#endif
  double CellShapingSpeed;

  int CPU_TimeBinCountMeasurements[TIMEBINS];
  double CPU_TimeBinMeasurements[TIMEBINS][NUMBER_OF_MEASUREMENTS_TO_RECORD];

  double MinGasHsmlFractional,	/**< minimum allowed SPH smoothing length in units of SPH gravitational
				   softening length */
    MinGasHsml;			/**< minimum allowed SPH smoothing length */

  double SofteningGas,		/**< for type 0 */
    SofteningHalo,		/**< for type 1 */
    SofteningDisk,		/**< for type 2 */
    SofteningBulge,		/**< for type 3 */
    SofteningStars,		/**< for type 4 */
    SofteningBndry;		/**< for type 5 */

  double SofteningGasMaxPhys,	/**< for type 0 */
    SofteningHaloMaxPhys,	/**< for type 1 */
    SofteningDiskMaxPhys,	/**< for type 2 */
    SofteningBulgeMaxPhys,	/**< for type 3 */
    SofteningStarsMaxPhys,	/**< for type 4 */
    SofteningBndryMaxPhys;	/**< for type 5 */

  double SofteningTable[6];	/**< current (comoving) gravitational softening lengths for each particle type */
  double ForceSoftening[6];	/**< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */

  double MassTable[6];

  /* some filenames */
  char InitCondFile[MAXLEN_PATH],
    OutputDir[MAXLEN_PATH],
    SnapshotFileBase[MAXLEN_PATH],
    EnergyFile[MAXLEN_PATH],
    CpuFile[MAXLEN_PATH],
    InfoFile[MAXLEN_PATH], TimingsFile[MAXLEN_PATH], RestartFile[MAXLEN_PATH], ResubmitCommand[MAXLEN_PATH], OutputListFilename[MAXLEN_PATH];

  /** table with desired output times */
  double OutputListTimes[MAXLEN_OUTPUTLIST];
  char OutputListFlag[MAXLEN_OUTPUTLIST];
  int OutputListLength;		/**< number of times stored in table of desired output times */


#ifdef USE_SFR			/* star formation and feedback sector */
  double CritOverDensity;
  double CritPhysDensity;
  double OverDensThresh;
  double PhysDensThresh;
  double TemperatureThresh;
  double EgySpecSN;
  double EgySpecCold;
  double FactorEVP;
  double TempSupernova;
  double TempClouds;
#endif

  double MaxSfrTimescale;
  double FactorSN; 

#if defined(REFINEMENT_SPLIT_CELLS) || defined(USE_SFR)
  MyIDType MaxID;
#endif
}
All;

extern struct particle_data
{
  MyDouble Pos[3] __attribute__((__aligned__(16)));
  MyDouble Mass;
  MyFloat  Vel[3] __attribute__((__aligned__(16)));
  MyIDType ID;

  //float GravCost[GRAVCOSTLEVELS];	/**< weight factors used for balancing the work-load */
  short int Type;		/**< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
  //short int TimeBin;
}
 *P,				/**< holds particle data on local processor */
 *DomainPartBuf;		/**< buffer for particle data used in domain decomposition */

extern struct sph_particle_data
{
  /* conserved variables */
  MyFloat Energy;
  MyFloat Momentum[3];
  MyFloat Volume;
  MyFloat OldMass;

  /* primitive variables */
  MyFloat Density;
  MyFloat Pressure;		/**< current pressure */
  MyFloat Utherm;

  /* variables for mesh  */
  MyFloat Center[3];		/* center of mass of cell */
  MyFloat VelVertex[3];		/* current vertex velocity (primitive variable) */
  MyFloat MaxDelaunayRadius;
  MyFloat Hsml;		        /* auxiliary search radius for points around a delaunay triangle */
  MyFloat SurfaceArea;
  MyFloat ActiveArea;

#ifdef COOLING
  MyFloat Ne;
#endif

#ifdef USE_SFR
  MyFloat Sfr;
#endif

  struct grad_data Grad;

#ifdef SECOND_DERIVATIVES
  struct hessian_data Hessian;
#endif

#ifdef VORONOI_DYNAMIC_UPDATE
  int first_connection;
  int last_connection;
#endif
}
 *SphP,				/**< holds SPH particle data on local processor */
 *DomainSphBuf;			/**< buffer for SPH particle data in domain decomposition */

extern struct NODE
{
  union
  {
    int suns[8];		/**< temporary pointers to daughter nodes */
    struct
    {
      MyFloat s[3] __attribute__((__aligned__(16)));		/**< center of mass of node */
      MyFloat mass;		/**< mass of node */
      float maxsoft;		/**< hold the maximum gravitational softening of particles */
      int sibling;		
      int nextnode;		
      /** The parent node of the node. (Is -1 for the root node.) */
      int father;		
    }
    d;
  }
  u;

  float center[3];           /**< geometrical center of node */
  float len;                 /**< sidelength of treenode */
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

extern int *Ngb_DomainNodeIndex;
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
      float range_min[3];
      int sibling;
      float range_max[3];
      int nextnode;
    }
    d;
  }
  u;
}
 *Ngb_Nodes; 
 
/*
 * proto.h derived (July 2012)
 */

extern "C" { 
void begrun1(void);
void endrun(int);
int init(void); 
void read_ic(const char *fname, int); 
 
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
}
 
#endif //AREPO_RT_ALLVARS_H
