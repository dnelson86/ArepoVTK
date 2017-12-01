/*
 * fileio_img.h
 * dnelson
 */

#ifndef AREPO_RT_FILEIO_IMG_H
#define AREPO_RT_FILEIO_IMG_H

#include <string>
using std::string;

void WriteImageTGA(const string &name, float *pixels, float *alpha, int xRes, int yRes,
                   int totalXRes, int totalYRes, int xOffset, int yOffset);
//void RGBSpectrum *ReadImageTGA(const string &name, int *w, int *h);

void WriteImagePNG(const string &name, float *pixels, float *alpha, int xRes, int yRes,
                   int totalXRes, int totalYRes, int xOffset, int yOffset);

#endif // AREPO_RT_FILEIO_IMG_H
