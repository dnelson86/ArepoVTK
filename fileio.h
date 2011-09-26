/*
 * fileio.h
 * dnelson
 */
 
#ifndef AREPO_RT_FILEIO_H
#define AREPO_RT_FILEIO_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
using namespace std;

#include <string>
using std::string;

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

    // data for public access
    string imageFile, rawRGBFile;
		string filename, paramFilename;
		
    int nTasks, nCores;
    bool quickRender, verbose, openWindow;
		
		int imageXPixels, imageYPixels;
		float swScale; // screenWindow mult factor * [-1,1]
		
		float cameraFOV;
		float cameraPosition[3], cameraLookAt[3];
		
		vector<string> tfSet;		
		
		bool drawBBox, drawTetra, drawVoronoi;
		bool projColDens, useDensGradients;		
		
		double viStepSize;
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

#include "ArepoRT.h"

#include <ctype.h>
#include <stdlib.h>

bool ReadFloatFile(const char *filename, vector<float> *values);
bool WriteFloatFile(const char *filename, float *values, int nx, int ny);
int parseSceneFile(const string &filename, int &nx, int &ny, int &nz, vector<float> *data);

void WriteImage(const string &name, float *pixels, float *alpha, int XRes, int YRes, 
                int totalXRes, int totalYRes, int xOffset, int yOffset);
								
#endif
