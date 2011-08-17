/*
 * fileio.h
 * dnelson
 */
 
#ifndef AREPO_RT_FILEIO_H
#define AREPO_RT_FILEIO_H

#include "ArepoRT.h"

#include <ctype.h>
#include <stdlib.h>

bool ReadFloatFile(const char *filename, vector<float> *values);
bool WriteFloatFile(const char *filename, float *values, int nx, int ny);
int parseSceneFile(const string &filename, int &nx, int &ny, int &nz, vector<float> *data);

void WriteImage(const string &name, float *pixels, float *alpha, int XRes, int YRes, 
                int totalXRes, int totalYRes, int xOffset, int yOffset);

#endif
