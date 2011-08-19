/*
 * ArepoRT.h
 * dnelson
 */
 
#ifndef AREPO_RT_H
#define AREPO_RT_H

// defines
#define AREPO_RT_VERSION   0.1
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

// standard libraries
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
class VolumeIntegrator;

// configuration options
struct ConfigStruct {
    ConfigStruct()
		{
		  nCores           = 0;
      quickRender      = quiet = openWindow = verbose = false;
      imageFile        = "frame.tga";
			rawRGBFile       = "frame.raw.txt";
			filename         = "test/Arepo2b.hdf5";
			paramFilename    = "test/param.txt";
			imageXPixels     = 1000; // 1024, 1920
			imageYPixels     = 1000; // 768, 1080
			viStepSize       = 0.1; // volume integration sub-stepping size (0=disabled)
			swScale          = 0.52f; // 0.52 for face on ortho with small border
			
			rgbEmit[0]    = 0.1f;   rgbEmit[1]    = 0.0f;   rgbEmit[2]    = 0.0f;
			rgbLine[0]    = 0.0f;   rgbLine[1]    = 0.0f;   rgbLine[2]    = 0.2f;
      rgbTetra[0]   = 0.1f;   rgbTetra[1]   = 0.1f;   rgbTetra[2]   = 0.1f;
      rgbVoronoi[0] = 0.00f;   rgbVoronoi[1] = 0.05f;   rgbVoronoi[2] = 0.00f;

		}
		
    int nCores;
		int imageXPixels, imageYPixels;
		double viStepSize;
		float rgbEmit[3], rgbLine[3], rgbTetra[3], rgbVoronoi[3];
		float swScale; // screenWindow mult factor * [-1,1]
    bool quickRender, quiet, verbose, openWindow;
    string imageFile, rawRGBFile;
		string filename, paramFilename;
};

extern struct ConfigStruct Config;

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

