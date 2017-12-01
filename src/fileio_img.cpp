/*
 * fileio_img.cpp
 * dnelson
 */
 
#include "fileio_img.h"
#include "ArepoRT.h"
#include "spectrum.h"

#define TO_BYTE(v) (uint8_t (Clamp(255.f   * powf((v), 1.0f/2.3f), 0.0f, 255.0f)))
#define TO_WORD(v) (uint16_t(Clamp(65535.f * powf((v), 1.0f/2.3f), 0.0f, 65535.0f)))

// TGA:
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <stdint.h>
typedef uint8_t byte;

typedef unsigned char uchar;

// PNG:
#include <png.h>

/**
 * gl_tga.c: TGA file format (TARGA) reader/writer.
 * License: GPL
 * Online License Link: http://www.gnu.org/licenses/gpl.html
 *
 * Copyright (c) 2003-2009 Jaakko Keranen <jaakko.keranen@iki.fi>
 * Copyright (c) 2009 Daniel Swanson <danij@dengine.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

// MACROS ------------------------------------------------------------------

#undef SHORT
#ifdef __BIG_ENDIAN__
#define SHORT(x)            shortSwap(x)
# else // Little-endian.
#define SHORT(x)            (x)
#endif

// TYPES -------------------------------------------------------------------

typedef struct {
    uchar           idLength; // Identification field size in bytes.
    uchar           colorMapType; // Type of the color map.
    uchar           imageType; // Image type code.
} tga_header_t;

// Color map specification.
typedef struct {
    int16_t         index; // Index of first color map entry.
    int16_t         length; // Number of color map entries.
    uchar           entrySize; // Number of bits in a color map entry (16/24/32).
} tga_colormapspec_t;

// Image specification.
typedef struct {
    int16_t         xOrigin; // X coordinate of lower left corner.
    int16_t         yOrigin; // Y coordinate of lower left corner.
    int16_t         width; // Width of the image in pixels.
    int16_t         height; // Height of the image in pixels.
    uchar           pixelDepth; // Number of bits in a pixel (16/24/32).
    uchar           attributeBits;
} tga_imagespec_t;


#ifdef __BIG_ENDIAN__
static int16_t shortSwap(int16_t n)
{
    return ((n & 0xff) << 8) | ((n & 0xff00) >> 8);
}
#endif

static bool writeByte(FILE* f, uchar b)
{
    return (fwrite(&b, 1, 1, f) == 1);
}

static bool writeShort(FILE* f, int16_t s)
{
    int16_t             v = SHORT(s);
    return (fwrite(&v, sizeof(v), 1, f) == 1);
}

/* static uchar readByte(FILE* f)
{
    uchar               v;
    fread(&v, sizeof(v), 1, f);
    return v;
}

static int16_t readShort(FILE* f)
{
    int16_t             v;
    fread(&v, sizeof(v), 1, f);
    return v;
} */

/**
 * @param idLength      Identification field size in bytes (max 255).
 *                      @c 0 indicates no identification field.
 * @param colorMapType  Type of the color map, @c 0 or @c 1:
 *                      @c 0 = color map data is not present.
 *                      @c 1 = color map data IS present.
 * @param imageType     Image data type code, one of:
 *                      @c 0 = no image data is present.
 *                      @c 1 = uncompressed, color mapped image.
 *                      @c 2 = uncompressed, true-color image.
 *                      @c 3 = uncompressed, grayscale image.
 *                      @c 9 = run-length encoded, color mapped image.
 *                      @c 10 = run-length encoded, true-color image.
 *                      @c 11 = run-length encoded, grayscale image.
 * @param file          Handle to the file to be written to.
 */
static void writeHeader(uchar idLength, uchar colorMapType, uchar imageType,
                        FILE* file)
{
    writeByte(file, idLength);
    writeByte(file, colorMapType? 1 : 0);
    writeByte(file, imageType);
}

/* static void readHeader(tga_header_t* dst, FILE* file)
{
    dst->idLength = readByte(file);
    dst->colorMapType = readByte(file);
    dst->imageType = readByte(file);
} */

/**
 * @param index         Index of first color map entry.
 * @param length        Total number of color map entries.
 * @param entrySize     Number of bits in a color map entry; 15/16/24/32.
 * @param file          Handle to the file to be written to.
 */
static void writeColorMapSpec(int16_t index, int16_t length,
                              uchar entrySize, FILE* file)
{
    writeShort(file, index);
    writeShort(file, length);
    writeByte(file, entrySize);
}

/* static void readColorMapSpec(tga_colormapspec_t* dst, FILE* file)
{
    dst->index = readShort(file);
    dst->length = readShort(file);
    dst->entrySize = readByte(file);
} */

/**
 * @param xOrigin       X coordinate of lower left corner.
 * @param yOrigin       Y coordinate of lower left corner.
 * @param width         Width of the image in pixels.
 * @param height        Height of the image in pixels.
 * @param pixDepth      Number of bits per pixel, one of; 16/24/32.
 * @param file          Handle to the file to be written to.
 */
static void writeImageSpec(int16_t xOrigin, int16_t yOrigin,
                           int16_t width, int16_t height, uchar pixDepth,
                           FILE* file)
{
    writeShort(file, xOrigin);
    writeShort(file, yOrigin);
    writeShort(file, width);
    writeShort(file, height);
    writeByte(file, pixDepth);

    /**
     * attributeBits:4; // Attribute bits associated with each pixel.
     * reserved:1; // A reserved bit; must be 0.
     * screenOrigin:1; // Location of screen origin; must be 0.
     * dataInterleave:2; // TGA_INTERLEAVE_*
     */
    writeByte(file, 0);
}

/* static void readImageSpec(tga_imagespec_t* dst, FILE* file)
{
    dst->xOrigin = readShort(file);
    dst->yOrigin = readShort(file);
    dst->width = readShort(file);
    dst->height = readShort(file);
    dst->pixelDepth = readByte(file);
    dst->attributeBits = readByte(file);
} */

/**
 * Save the rgb8888 buffer as Targa 24.
 *
 * @param filename      Path to the file to be written to (need not exist).
 * @param w             Width of the image in pixels.
 * @param h             Height of the image in pixels.
 * @param buf           Ptr to the image data to be written.
 *
 * @return              Non-zero iff successful.
 */
void WriteImageTGA(const string &name, float *pixels,
        float *alpha, int xRes, int yRes,
        int totalXRes, int totalYRes,
        int xOffset, int yOffset)
{
    FILE*               file;
    uchar*              outBuf;

    if ((file = fopen(name.c_str(), "wb")) == NULL) {
        cout << "WriteImageTGA: Unable to open output file!" << endl;
        return;
    }

    // No identification field, no color map, Targa type 2 (unmapped RGB).
    writeHeader(0, 0, 2, file);
    writeColorMapSpec(0, 0, 0, file);
    writeImageSpec(0, 0, xRes, yRes, 24, file);

    // The save format is BGR.
    outBuf = (uchar *)malloc(xRes * yRes * 3);
    uchar *dst = outBuf;
    for (int y = yRes-1; y >= 0; --y) {
        for (int x = 0; x < xRes; ++x) {
            dst[0] = TO_BYTE(pixels[3*(y*xRes+x)+2]);
            dst[1] = TO_BYTE(pixels[3*(y*xRes+x)+1]);
            dst[2] = TO_BYTE(pixels[3*(y*xRes+x)+0]);
            dst += 3;
        }
    }
    if (fwrite(outBuf, 1, 3 * xRes * yRes, file) != uint32_t(3*xRes*yRes))
        cout << "WriteImageTGA: Error writing file!" << endl;
    free(outBuf);
    fclose(file);
		
		if( Config.writeRGBA8bit || Config.writeRGB16bit ) {
			cout << "WARNING: Requested RGBA 8bit or RGB 16bit, neither implemented yet for TGA!" << endl;
		}
}

/**
 * Loads a TGA image (not the RLE types though)
 */
 /*
void RGBSpectrum *ReadImageTGA(const string &name, int *width, int *height)
{
    int                 x, y, pixbytes;
    tga_header_t        header;
    tga_colormapspec_t  colorMapSpec;
    tga_imagespec_t     imageSpec;
    uchar*              srcBuf;
    const uchar*        src;

    FILE *file = fopen(name.c_str(), "rb");
    if (!file) {
        cout << "ReadImageTGA: Unable to open input file!" << endl;
        return NULL;
    }

    // Read and check the header.
    readHeader(&header, file);
    readColorMapSpec(&colorMapSpec, file);
    readImageSpec(&imageSpec, file);

    if (((imageSpec.attributeBits & 0xf) != 8 &&  // num attribute bits
         (imageSpec.attributeBits & 0xf) != 0) ||
        ((imageSpec.attributeBits & 0xc0) != 0) || // no interleaving
        (header.imageType == 2 &&
          (imageSpec.pixelDepth != 32 && imageSpec.pixelDepth != 24)) ||
        (header.imageType == 3 &&
          (imageSpec.pixelDepth != 8)) ||
        (header.imageType != 2 && header.imageType != 3)) {
        cout << "ReadImageTGA: Unrecognized format!" << "type=" <, header.imageType 
				     << " pxsize=" << imageSpec.pixelDepth << " abits=" << imageSpec.attributeBits << ")" << endl;
        fclose(file);
        return NULL;
    }

    *width = imageSpec.width;
    *height = imageSpec.height;

    // Determine format.
    if (imageSpec.pixelDepth == 32)
        pixbytes = 4;
    else if (imageSpec.pixelDepth == 24)
        pixbytes = 3;
    else if (imageSpec.pixelDepth == 8) {
        pixbytes = 1;
    } else {
				cout << "ReadImageTGA: pixelDepth error!" << endl;
				return NULL;
		}

    // Read the pixel data.
    int size = *width * *height * pixbytes;
    srcBuf = (uchar *)malloc(size);
    if (fread(srcBuf, 1, size, file) != (uint32_t)size) {
        cout << "ReadImageTGA: Premature end of file!" << endl;
        free(srcBuf);
        fclose(file);
        return NULL;
    }

    // "Unpack" the pixels (origin in the lower left corner).
    // TGA pixels are in BGRA format.
    src = srcBuf;
    RGBSpectrum *ret = new RGBSpectrum[*width * *height];
    RGBSpectrum *dst = ret;
    for (y = *height - 1; y >= 0; y--)
        for (x = 0; x < *width; x++) {
            if (pixbytes == 1)
                *dst++ = RGBSpectrum((*src++) / 255.f);
            else {
                float c[3];
                c[2] = (*src++) / 255.f;
                c[1] = (*src++) / 255.f;
                c[0] = (*src++) / 255.f;
                *dst++ = RGBSpectrum::FromRGB(c);
                if (pixbytes == 4)
                    ++src;
            }
    }

    bool flipH = ((imageSpec.attributeBits & 0x10) == 0x10);
    bool flipV = ((imageSpec.attributeBits & 0x20) == 0x20);
    if (flipH) {
        for (y = 0; y < *height; ++y)
            for (x = 0; x < *width / 2; ++x)
                swap(ret[y * *width + x], ret[y * *width + (*width - 1 - x)]);
    }
    if (flipV) {
        for (y = 0; y < *height/2; ++y)
            for (x = 0; x < *width; ++x)
                swap(ret[y * *width + x], ret[(*height - 1 - y) * *width + x]);
    }
    free(srcBuf);
    fclose(file);
    cout << "ReadTGAImage: Read " << name << " (" << *width << " x " << *height << ")" << endl;
    return ret;
}
*/

void pbrt_png_error(png_structp png_, png_const_charp msg)
{
   cout << "libpng error: " << msg << endl;
}

void WriteImagePNG(const string &name, float *pixels, float *alpha, int xRes, int yRes,
                   int totalXRes, int totalYRes, int xOffset, int yOffset)
{
	// PNG file
	unsigned char** rows;
	png_structp png;
	png_infop info;
	FILE *fp;
	
	png_text text;
	text.compression = PNG_TEXT_COMPRESSION_NONE;
	text.key = (png_charp) "Software";
	text.text = (png_charp) "ArepoVTK";
	text.text_length = 8;
		
	// write 8-bit RGB
	// ---------------
	if( Config.writeRGB8bit )
	{
		// open file
		if ((fp = fopen(name.c_str(), "wb")) == NULL) {
				cout << "WriteImagePNG: Unable to open output file!" << endl;
				return;
		}
		
		// write PNG header
		png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, (png_error_ptr) pbrt_png_error, NULL);
		info = png_create_info_struct(png);
		png_init_io(png, fp);
		png_set_text(png, info, &text, 1);

		//png_color_16 black = {0};
		//png_set_background(png, &black, PNG_BACKGROUND_GAMMA_SCREEN, 0, 255.0);
		
		// gAMA: (unlike TGA we do not explicitly gamma scale the discretized values before writing)
		//float screenGamma = 1.0;
		//float fileGamma   = 1.0;
		//png_set_gamma(png, screenGamma, fileGamma);
		
		// cHRM: is: CIE x,y chromacities of R, G, B and white.
		//
		// x = X / (X + Y + Z)
		// y = Y / (X + Y + Z)
		/*
		float rgbWhite[3] = {1, 1, 1};
		float rgbRed[3] = {1, 0, 0};
		float rgbGreen[3] = {0, 1, 0};
		float rgbBlue[3] = {0, 0, 1};
		float xyzWhite[3];
		float xyzRed[3];
		float xyzGreen[3];
		float xyzBlue[3];

		Spectrum::FromRGB(rgbWhite).ToXYZ(xyzWhite);
		Spectrum::FromRGB(rgbRed).ToXYZ(xyzRed);
		Spectrum::FromRGB(rgbGreen).ToXYZ(xyzGreen);
		Spectrum::FromRGB(rgbBlue).ToXYZ(xyzBlue);
		
		float whiteX = xyzWhite[0] / (xyzWhite[0] + xyzWhite[1] + xyzWhite[2]);
		float whiteY = xyzWhite[1] / (xyzWhite[0] + xyzWhite[1] + xyzWhite[2]);
		float redX = xyzRed[0] / (xyzRed[0] + xyzRed[1] + xyzRed[2]);
		float redY = xyzRed[1] / (xyzRed[0] + xyzRed[1] + xyzRed[2]);
		float greenX = xyzGreen[0] / (xyzGreen[0] + xyzGreen[1] + xyzGreen[2]);
		float greenY = xyzGreen[1] / (xyzGreen[0] + xyzGreen[1] + xyzGreen[2]);
		float blueX = xyzBlue[0] / (xyzBlue[0] + xyzBlue[1] + xyzBlue[2]);
		float blueY = xyzBlue[1] / (xyzBlue[0] + xyzBlue[1] + xyzBlue[2]);

		png_set_cHRM(png, info, whiteX, whiteY, redX, redY, greenX, greenY, blueX, blueY);	
		*/
		// IHDR
		png_set_IHDR(png, info, xRes, yRes, 8,
								 PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
								 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
					
		// do write
		rows = (png_bytep*) malloc(yRes * sizeof(png_bytep));
		rows[0] = (png_bytep) malloc(xRes * yRes * 3);
		for (int i = 1; i < yRes; i++)
				rows[i] = rows[0] + i * xRes * 3;

		for (int x = xOffset; x < xOffset + xRes; ++x)
		{
			for (int y = yOffset; y < yOffset + yRes; ++y)
			{
				char r = TO_BYTE( pixels[(x + y * totalXRes) * 3 + 0] );
				char g = TO_BYTE( pixels[(x + y * totalXRes) * 3 + 1] );
				char b = TO_BYTE( pixels[(x + y * totalXRes) * 3 + 2] );
				rows[y - yOffset][x * 3 + 0] = r;
				rows[y - yOffset][x * 3 + 1] = g;
				rows[y - yOffset][x * 3 + 2] = b;
			}
		}

		png_set_rows(png, info, rows);
		png_write_png(png, info, PNG_TRANSFORM_IDENTITY, NULL);
		
		fclose(fp);
		delete rows[0];
	}
	
	// write 8-bit RGBA
	// ----------------
	if( Config.writeRGBA8bit )
	{
		string name2 = name.substr(0,name.find(".png")) + "_rgba.png";
		if ((fp = fopen(name2.c_str(), "wb")) == NULL) {
				cout << "WriteImagePNG: Unable to open RGBA output file!" << endl;
				return;
		}
		
		// PNG header
		png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, (png_error_ptr) pbrt_png_error, NULL);
		info = png_create_info_struct(png);
		png_init_io(png, fp);
		png_set_text(png, info, &text, 1);

		// IHDR
		png_set_IHDR(png, info, xRes, yRes, 8,
								 PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
								 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
					
		// do write
		rows = (png_bytep*) malloc(yRes * sizeof(png_bytep));
		rows[0] = (png_bytep) malloc(xRes * yRes * 4);
		for (int i = 1; i < yRes; i++)
				rows[i] = rows[0] + i * xRes * 4;

		for (int x = xOffset; x < xOffset + xRes; ++x)
		{
			for (int y = yOffset; y < yOffset + yRes; ++y)
			{
				char r = TO_BYTE( pixels[(x + y * totalXRes) * 3 + 0] );
				char g = TO_BYTE( pixels[(x + y * totalXRes) * 3 + 1] );
				char b = TO_BYTE( pixels[(x + y * totalXRes) * 3 + 2] );
				char a = TO_BYTE( 1.0 - alpha[x + y * totalXRes] );
				//char a = TO_BYTE( 1.0 );
				rows[y - yOffset][x * 4 + 0] = r;
				rows[y - yOffset][x * 4 + 1] = g;
				rows[y - yOffset][x * 4 + 2] = b;
				rows[y - yOffset][x * 4 + 3] = a;
			}
		}

		png_set_rows(png, info, rows);
		png_write_png(png, info, PNG_TRANSFORM_IDENTITY, NULL);
		
		fclose(fp);
		delete rows[0];
	}
	
	// write 16-bit RGB
	// ----------------
	if( Config.writeRGB16bit )
	{
		string name3 = name.substr(0,name.find(".png")) + "_16bit.png";
		if ((fp = fopen(name3.c_str(), "wb")) == NULL) {
				cout << "WriteImagePNG: Unable to open 16bit output file!" << endl;
				return;
		}
		
		// PNG header
		png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, (png_error_ptr) pbrt_png_error, NULL);
		png_set_swap(png);
		info = png_create_info_struct(png);
		png_init_io(png, fp);
		png_set_text(png, info, &text, 1);

		// IHDR
		png_set_IHDR(png, info, xRes, yRes, 16,
								 PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
								 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
					
		// do write
		rows = (png_bytep*) malloc(yRes * sizeof(png_bytep));
		rows[0] = (png_bytep) malloc(xRes * yRes * 3 * 2);
		for (int i = 1; i < yRes; i++)
				rows[i] = rows[0] + i * xRes * 3 * 2;

		for (int x = xOffset; x < xOffset + xRes; ++x)
		{
			for (int y = yOffset; y < yOffset + yRes; ++y)
			{
				uint16_t r = TO_WORD( pixels[(x + y * totalXRes) * 3 + 0] );
				uint16_t g = TO_WORD( pixels[(x + y * totalXRes) * 3 + 1] );
				uint16_t b = TO_WORD( pixels[(x + y * totalXRes) * 3 + 2] );
				//uint16_t a = TO_WORD( 1.0 - alpha[x + y * totalXRes] );
				//uint16_t a = TO_WORD( 1.0 );
				// TODO: alpha has block structure corresponding to tasks, needs fix (in normalization?)
				
				// swap endianness and write
				rows[y - yOffset][x * 3 * 2 + 0] = ((uint8_t*) &r)[1];
				rows[y - yOffset][x * 3 * 2 + 1] = ((uint8_t*) &r)[0];
				rows[y - yOffset][x * 3 * 2 + 2] = ((uint8_t*) &g)[1];
				rows[y - yOffset][x * 3 * 2 + 3] = ((uint8_t*) &g)[0];
				rows[y - yOffset][x * 3 * 2 + 4] = ((uint8_t*) &b)[1];
				rows[y - yOffset][x * 3 * 2 + 5] = ((uint8_t*) &b)[0];
			}
		}

		png_set_rows(png, info, rows);
		png_write_png(png, info, PNG_TRANSFORM_IDENTITY, NULL);
		
		fclose(fp);
		delete rows[0];
	}
	
}