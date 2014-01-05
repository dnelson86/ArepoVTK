/*
 * ArepoRT.h
 * dnelson
 */
 
#ifndef AREPO_RT_H
#define AREPO_RT_H

// defines
#define AREPO_RT_VERSION    0.41
#define L1_CACHE_LINE_SIZE  64
#define FILTER_TABLE_SIZE   16
#define TASK_MULT_FACT      8 //32
#define TASK_MAX_PIXEL_SIZE 100 //16
#define INFINITY            FLT_MAX
#define INSIDE_EPS          1.0e-6
#define AUXMESH_ALLOC_SIZE  1000

#define MSUN_PER_PC3_IN_CGS 6.769e-23

// for selective load (temporary)
#define READ_PARTTYPE 0
#define TILESIZE      256
#define FILL_NUMPART  1000

// behavior options

//#define USE_LINEALGO_BRESENHAM
#define USE_LINEALGO_WU

/* interpolation methods (choose one) */

//#define NATURAL_NEIGHBOR_INTERP
//#define NATURAL_NEIGHBOR_IDW
#define NATURAL_NEIGHBOR_SPHKERNEL
//#define NNI_WATSON_SAMBRIDGE
//#define NNI_LIANG_HALE
//#define DTFE_INTERP
//#define CELL_GRADIENTS_ONLY

/* interpolation method options */

//#define NO_GHOST_CONTRIBS // only for SPHKERNEL, do not use
                            // ghosts for hsml/TF (i.e. for reflective BCs but we are doing
                            // periodic meshing, the ghosts are incorrect and should be skipped)
//#define NNI_DISABLE_EXACT // for bruteforce NNI disable exact geometry computations
//#define NATURAL_NEIGHBOR_INNER // for IDW or SPHKERNEL do neighbors of neighbors
//#define BRUTE_FORCE            // for IDW or SPHKERNEL, calculate over all NumGas in the box

/* exponent of distance, greater values assign more influence to values closest to the 
 * interpolating point, approaching piecewise constant for large POWER_PARAM.
 * in N dimensions, if p <= N, the interpolated values are dominated by points far away,
 * which is rather bizarre. note: p is actually 2p since we skip the sqrt.
 */
#define POWER_PARAM 2.0

/* use <1 for more smoothing, makes hsml_used bigger than hsml_ngb_max
 * use >1 for less smoothing, make hsml_used smaller than hsml_ngb_max 
 *   (at least some some neighbors will not contribute, and if this is too big, the
 *    containing cell may also not contribute, leading to black holes)
 */
#define HSML_FAC 1.0

/* special behavior */
//#define DEBUG_VERIFY_INCELL_EACH_STEP
//#define DEBUG_VERIFY_ENTRY_CELLS
//#define DISABLE_MEMORY_MANAGER

//#define DUMP_VORONOI_MESH
//#define DUMP_MESH_TEXT

#ifdef DEBUG
#define IF_DEBUG(expr) (expr)
#else
#define IF_DEBUG(expr) ((void)0)
#endif

// includes
#include <math.h>
#include <sys/time.h>
#include <stdint.h>
#include <malloc.h> //memalign

#include <typeinfo>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
using namespace std;

#include <string>
using std::string;
#include <vector>
using std::vector;
#include <map>

#include "fileio.h" //Config

class ConfigSet;
extern ConfigSet Config;

// global forward declarations
template <typename T, int logBlockSize = 2> class BlockedArray;
struct Matrix4x4;

class Vector;
class Point;
class Ray;
class BBox;

class Renderer;
class Transform;
struct Intersection;

template <int nSamples> class CoefficientSpectrum;
class RGBSpectrum;
typedef RGBSpectrum Spectrum;

class Camera;
class Film;
class RNG;
class Sampler;
struct CameraSample;
struct Sample;

class Arepo;
class ArepoMesh;
class ArepoTree;
class ArepoSnapshot;
class Scene;
class VolumeRegion;
class DensityRegion;
class VolumeGridDensity;
class TransferFunction;
class VolumeIntegrator;
class VoronoiIntegrator;
class TreeSearchIntegrator;
class Timer;

// util functions
template <typename T> string toStr(T num)
{
	stringstream ss;
	ss << num;
	return ss.str();
}

// global inlines
inline float Lerp(float t, float v1, float v2) {
    return (1.0f - t) * v1 + t * v2;
}

inline float Radians(float deg) {
    return ((float)M_PI/180.0f) * deg;
}
inline float Degrees(float rad) {
    return (180.0f/(float)M_PI) * rad;
}

/*inline float Clamp(float val, float low, float high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}*/
inline double Clamp(double val, double low, double high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}
inline int Clamp(int val, int low, int high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}
inline unsigned int RoundUpPowerOfTwo(int val) {
    val--;
    val |= val >> 1; val |= val >> 2; val |= val >> 4; val |= val >> 8; val |= val >> 16;
    val++;
    return val;
}

#endif

