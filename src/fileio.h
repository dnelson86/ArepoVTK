/*
 * fileio.h
 * dnelson
 */
 
#ifndef AREPO_RT_FILEIO_H
#define AREPO_RT_FILEIO_H

#include "ArepoRT.h"

// configuration options
class ConfigSet {
public:
  // construction
  ConfigSet() {
      delim   = "=";
      comment = "%";
  }

  // methods
  void ReadFile(string cfgfile);
  void print();
  
  template<class T> T readValue(const string &key) const;
  template<class T> T readValue(const string &key, const T &defaultValue) const;
  
  void splitStrArray(const string str, float *rgb); // size 3

  // Input/Output
  string imageFile, rawRGBFile;
  string filename, paramFilename;
  bool writeRGB8bit, writeRGB16bit, writeRGBA8bit;

  bool dumpMeshText, dumpMeshBinary, dumpMeshCells;
  
  // General
  int nTasks, nCores;
  bool quickRender, verbose, openWindow;
  
  int totNumJobs, curJobNum;
  int jobExpansionFac, expandedJobNum;
  string maskFileBase;
  float maskPadFac;
  
  // Frame/Camera
  int imageXPixels, imageYPixels;
  float swScale; // screenWindow mult factor * [-1,1]
  
  string cameraType;
  float cameraFOV;
  float cameraPosition[3], cameraLookAt[3], cameraUp[3];
  
  // Data Processing
  float recenterBoxCoords[3];
  bool convertUthermToKelvin;
  bool takeLogUtherm, takeLogDens;
  
  // Transfer Functions
  int readPartType;
  vector<string> tfSet;
  
  // Animation
  vector<string> kfSet; // key frames
  
  int startFrame, numFrames;
  float timePerFrame;
  float minScale, maxScale; // scale all RGB values from [min,max]->[0,1]
  float minAlpha, maxAlpha; // scale raw density integrals from [min,max]->[0,1] and write as alpha channel
  
  // Render
  bool drawBBox, drawTetra, drawVoronoi, drawSphere;
  bool projColDens;   
  
  int nTreeNGB;
  float viStepSize;
  float rayMaxT;
  float rgbLine[3], rgbTetra[3], rgbVoronoi[3];
  float rgbAbsorb[3];

private:
  // for reading config file
  map<string,string> parsedParams;
  
  string delim, comment;
  
  typedef map<string,string>::iterator       map_i;
  typedef map<string,string>::const_iterator map_ci;
  
  template<class T> static T string_to_T(const string& str);
};

bool ReadFloatFile(const char *filename, vector<float> *values);
bool WriteFloatFile(const char *filename, float *values, int nx, int ny);
int parseSceneFile(const string &filename, int &nx, int &ny, int &nz, vector<float> *data);

int loadDiscreteColorTable(const string &filename, vector<float> *colorTableVals);

void WriteImage(const string &name, float *pixels, float *alpha, int XRes, int YRes, 
                int totalXRes, int totalYRes, int xOffset, int yOffset);
                
#endif
