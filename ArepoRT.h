/*
 * ArepoRT.h
 * dnelson
 */
 
#ifndef AREPO_RT_H
#define AREPO_RT_H

// defines
#define AREPO_RT_VERSION   0.2
#define L1_CACHE_LINE_SIZE 64
#define FILTER_TABLE_SIZE  16
#define INFINITY           FLT_MAX
#define INSIDE_EPS         1.0e-6

//#define USE_LINEALGO_BRESENHAM
#define USE_LINEALGO_WU

#ifdef DEBUG
#define IF_DEBUG(expr) (expr)
#else
#define IF_DEBUG(expr) ((void)0)
#endif

#ifdef ENABLE_AREPO
#define IF_AREPO(expr) (expr)
#else
#define IF_AREPO(expr) ((void)0)
typedef int ArepoMesh; //pointers
#endif

// includes
#include <math.h>
#include <sys/time.h>

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include <string>
using std::string;
#include <vector>
using std::vector;

#include <malloc.h> //memalign

#include "fileio.h" //Config

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

#ifdef ENABLE_AREPO
class Arepo;
class ArepoMesh;
#endif

class Scene;
class VolumeRegion;
class DensityRegion;
class VolumeGridDensity;
class TransferFunction;

class VolumeIntegrator;
class VoronoiIntegrator;

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

inline float Clamp(float val, float low, float high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}
inline int Clamp(int val, int low, int high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}

#endif

