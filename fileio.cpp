/*
 * fileio.cpp
 * dnelson
 */

#include "fileio.h"

// ConfigSet
void ConfigSet::ReadFile(string cfgfile)
{
		cout << "Reading Configuration File [" << cfgfile << "]." << endl;
		
		ifstream is(cfgfile.c_str());
		if (!is) {
				cout << "ConfigSet: ERROR opening configuration file [" << cfgfile << "]." << endl;
				exit(1113);
		}
		
		// prepare to parse
		const string::size_type skip = delim.length();
		string nextLine = "";
		
		static const char whitespace[] = " \n\t\v\r\f";
		
		while (is || nextLine.length() > 0) {
				// read line
				string line = "";
				bool term   = false;
				
				if (nextLine.length() > 0) {
						line     = nextLine;
						nextLine = "";
				} else {
						getline(is,line);
				}
				
				// ignore inline and full line comments
				line = line.substr(0,line.find(comment));
				
				// split into key = value pair
				string::size_type offset = line.find(delim);
				
				if (offset < string::npos) {
						// get key
						string key = line.substr(0,offset);
						line.replace(0,offset+skip,"");
				
						// trim
						key.erase( 0, key.find_first_not_of(whitespace) );
						key.erase( key.find_last_not_of(whitespace) + 1U );
						
						line.erase( 0, line.find_first_not_of(whitespace) );
						line.erase( line.find_last_not_of(whitespace) + 1U );
				
						// store key,value (map type, keys unique)
						parsedParams[key] = line;
						//IF_DEBUG(cout << " [" << key << "] = " << line << endl);
				}
		}
		
		// fill public data members
		imageFile     = readValue<string>("imageFile",  "frame.tga");
		rawRGBFile    = readValue<string>("rawRGBFile", "frame.raw.txt");
		filename      = readValue<string>("filename");
		paramFilename = readValue<string>("paramFilename");
		
		nTasks        = readValue<int> ("nTasks",        1);
		quickRender   = readValue<bool>("quickRender",   false);
		openWindow    = readValue<bool>("openWindow",    false);
		verbose       = readValue<bool>("verbose",       false);
		
		imageXPixels  = readValue<int>  ("imageXPixels", 500);
		imageYPixels  = readValue<int>  ("imageYPixels", 500);
		swScale       = readValue<float>("swScale",      1.0f);
		
		drawBBox      = readValue<bool>("drawBBox",      true);
		drawTetra     = readValue<bool>("drawTetra",     true);	
		drawVoronoi   = readValue<bool>("drawVoronoi",   false);

		projColDens   = readValue<bool>("projColDens",   false);
		
		viStepSize    = readValue<float>("viStepSize",   0.0f); // disabled by default
		
		//TODO: temp rgb triplets input		
		splitStrArray( readValue<string>("rgbEmit",     "0.1  0.0  0.0")  , &rgbEmit[0]    );
		splitStrArray( readValue<string>("rgbLine",     "0.1  0.1  0.1")  , &rgbLine[0]    );
		splitStrArray( readValue<string>("rgbTetra",    "0.01 0.01 0.01") , &rgbTetra[0]   );
		splitStrArray( readValue<string>("rgbVoronoi",  "0.0  0.05 0.0")  , &rgbVoronoi[0] );

		splitStrArray( readValue<string>("rgbAbsorb",   "0.0  0.05 0.0")  , &rgbAbsorb[0]  );
		
		// basic validation
		if (projColDens && viStepSize) {
				cout << "Config: ERROR! projColDens and viStepSize are incompatible options!" << endl;
				exit(1120);
		}
}

void ConfigSet::print()
{
		map_i pi;
		
		// TODO: fix output of rgb triplets / store in some other way
		cout << endl << "CONFIGURATION PARAMETERS USED:" << endl << endl;
		for (pi = parsedParams.begin(); pi != parsedParams.end(); ++pi) {
				cout << " " << pi->first << " " << delim << " " << pi->second << endl;
		}
		cout << endl;

}

template<class T> T ConfigSet::readValue(const string &key) const
{
		// no default value = required key, die if missing
		map_ci pi = parsedParams.find(key);
		
		if( pi == parsedParams.end() ) {
				cout << " ERROR: Required parameter [" << key << "] missing from configuration file." << endl;
				exit(1114);
		}
		
		return string_to_T<T>( pi->second );
}

template<class T> T ConfigSet::readValue(const string &key, const T &defaultValue) const
{
		// default value supplied = optional param in config file
		map_ci pi = parsedParams.find(key);
		
		if( pi == parsedParams.end() )
				return defaultValue;
				
		return string_to_T<T>( pi->second );
}

template<class T> T ConfigSet::string_to_T(const string &str)
{
		T t;
		istringstream ist(str);
		
		ist >> t; // use >> operator to typecast to T
		return t;
}

template<> inline bool ConfigSet::string_to_T<bool>(const string &str)
{
		// special case conversion to bool
		bool ret   = true;
		string strUpper = str;
		
		// uppercase
		for(string::iterator p = strUpper.begin(); p != strUpper.end(); ++p)
				*p = toupper(*p);
			
		if(strUpper == string("FALSE") || strUpper == string("NO") || strUpper == string("0"))
				ret = false;
			
		return ret;
}

template<> inline string ConfigSet::string_to_T<string>(const string &str)
{
	// special case conversion to string, otherwise spaces get lost with >>
	return str;
}

void ConfigSet::splitStrArray(const string &str, float *rgb) //size=3
{
		char *pch = strtok((char*)str.c_str(), " ,");
		for (int i=0; i < 3; i++) {
				rgb[i] = atof(pch);
				pch = strtok(NULL," ,");
		}
}

// Image Output

static void WriteImageTGA(const string &name, float *pixels, float *alpha, int xRes, int yRes,
                          int totalXRes, int totalYRes, int xOffset, int yOffset);
//static RGBSpectrum *ReadImageTGA(const string &name, int *w, int *h);

bool ReadFloatFile(const char *filename, vector<float> *values)
{
    FILE *f = fopen(filename, "r");
    if (!f) {
        return false;
    }

    int c;
    bool inNumber = false;
    char curNumber[32];
    int curNumberPos = 0;
    int lineNumber = 1;
    while ((c = getc(f)) != EOF) {
        if (c == '\n') ++lineNumber;
        if (inNumber) {
            if (isdigit(c) || c == '.' || c == 'e' || c == '-' || c == '+')
                curNumber[curNumberPos++] = c;
            else {
                curNumber[curNumberPos++] = '\0';
                values->push_back(atof(curNumber));
                //Assert(curNumberPos < (int)sizeof(curNumber));
                inNumber = false;
                curNumberPos = 0;
            }
        }
        else {
            if (isdigit(c) || c == '.' || c == '-' || c == '+') {
                inNumber = true;
                curNumber[curNumberPos++] = c;
            }
            else if (c == '#') {
                while ((c = getc(f)) != '\n' && c != EOF)
                    ;
                ++lineNumber;
            }
            else if (!isspace(c)) {
								cout << "Unexpected text found at line " << lineNumber << " of float file "
								     << filename << endl;
            }
        }
    }
    fclose(f);
    return true;
}

bool WriteFloatFile(const char *filename, float *values, int nx, int ny)
{
		// c - text
    FILE *f = fopen(filename, "w");
    if (!f) {
        cout << "WriteFloatFile: Unable to open file " << filename << endl;
        return false;
		}
		
		IF_DEBUG(cout << "WriteFloatFile(" << filename << ") nx = " << nx << " ny = " << ny << endl);
		
		int y, x;
		int offset = 0;
		float amp = 0.0f;
		
		fprintf(f, "%d\n%d\n", nx, ny);
		
    for (y = 0; y < ny; ++y) {
        for (x = 0; x < nx; ++x) {
						fprintf(f, "(%f %f %f)", values[3*offset  ],
						                         values[3*offset+1],
																		 values[3*offset+2]);
						offset++;
				}
				fprintf(f, "\n");
		}
		
		fprintf(f, "amp\n");
		
		offset = 0;
		
    for (y = 0; y < ny; ++y) {
        for (x = 0; x < nx; ++x) {
						amp = sqrt( values[3*offset]*values[3*offset] + values[3*offset+1]*values[3*offset+1]
											+ values[3*offset+2]*values[3*offset+2] );
						fprintf(f, "%f ", amp);
						offset++;
				}
				fprintf(f, "\n");
		}
		
    fclose(f);
		
		// c++ - binary
/*		string fname = filename;
		fname += ".bin";
		ofstream f2(fname.c_str(), ofstream::binary);
		
		if (!f2.is_open()) {
        cout << "WriteFloatFile: Unable to open file " << filename << endl;
        return false;
		}
		
		cout << "WriteFloatFile(" << fname << ") nx = " << nx << " ny = " << ny << endl;		
		
		f2.write((const char*)&nx,sizeof(int));
		f2.write((const char*)&ny,sizeof(int));
		f2.write((const char*)&values,sizeof(float)*3*nx*ny);
		f2.close(); */

    return true;
}

int parseSceneFile(const string &filename, int &nx, int &ny, int &nz, vector<float> *data)
{
		int ni = 0;
		vector<float> vals;
		
		if (!ReadFloatFile(filename.c_str(), &vals)) {
				return 0;
		}

		// size
		if (vals.size() >= 4) {
			ni = (int)vals[0];
			nx = (int)vals[1];
			ny = (int)vals[2];
			nz = (int)vals[3];
		}

		// errorcheck size
		if (ni != nx*ny*nz) {
			cout << "Error: parseSceneFile ni = " << ni << " but nx*ny*nz = " << nx*ny*nz << endl;
			return 0;
		}
		if (ni != (int)ni || nx != (int)nx || ny != (int)ny || nz != (int)nz) {
			cout << "Error: parseSceneFile non-integer size, ni nx ny nz = " 
					 << ni << " " << nx << " " << ny << " " << nz << endl;
			return 0;
		}
	
#if defined(DEBUG) && 0
		cout << "parseSceneFile volumeGrid: ni = " << ni 
		              << " nx = " << nx << " ny = " << ny << " nz = " << nz << endl;
									
		cout << "parseSceneFile volumeGrid data:";
		for (int i=4; i < vals.size(); i++)
		  cout << " " << vals[i];
		cout << endl;
#endif
		
		for (int i=4; i < vals.size(); i++)
			data->push_back(vals[i]);
		
		return ni;
}

void WriteImage(const string &name, float *pixels, float *alpha, int xRes,
                int yRes, int totalXRes, int totalYRes, int xOffset, int yOffset)
{
    if (name.size() >= 5) {
        uint32_t suffixOffset = name.size() - 4;
#ifdef HAVE_OPENEXR
        if (!strcmp(name.c_str() + suffixOffset, ".exr") ||
            !strcmp(name.c_str() + suffixOffset, ".EXR")) {
             WriteImageEXR(name, pixels, alpha, xRes, yRes, totalXRes,
                           totalYRes, xOffset, yOffset);
             return;
        }
#endif
        if (!strcmp(name.c_str() + suffixOffset, ".tga") ||
            !strcmp(name.c_str() + suffixOffset, ".TGA")) {
            WriteImageTGA(name, pixels, alpha, xRes, yRes, totalXRes,
                          totalYRes, xOffset, yOffset);
            return;
        }
    }
    cout << "WriteImage: Can't determine image file type from suffix of filename '" << name << "'!" << endl;
}

// TGA Function Definitions
/**\file
 *\section License
 * License: GPL
 * Online License Link: http://www.gnu.org/licenses/gpl.html
 *
 *\author Copyright (c) 2003-2009 Jaakko Keranen <jaakko.keranen@iki.fi>
 *\author Copyright (c) 2009 Daniel Swanson <danij@dengine.net>
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

/**
 * gl_tga.c: TGA file format (TARGA) reader/writer.
 */

// HEADER FILES ------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <stdint.h>
typedef uint8_t byte;

typedef unsigned char uchar;

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

static uchar readByte(FILE* f)
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
}

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

static void readHeader(tga_header_t* dst, FILE* file)
{
    dst->idLength = readByte(file);
    dst->colorMapType = readByte(file);
    dst->imageType = readByte(file);
}

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

static void readColorMapSpec(tga_colormapspec_t* dst, FILE* file)
{
    dst->index = readShort(file);
    dst->length = readShort(file);
    dst->entrySize = readByte(file);
}

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

static void readImageSpec(tga_imagespec_t* dst, FILE* file)
{
    dst->xOrigin = readShort(file);
    dst->yOrigin = readShort(file);
    dst->width = readShort(file);
    dst->height = readShort(file);
    dst->pixelDepth = readByte(file);
    dst->attributeBits = readByte(file);
}

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
#define TO_BYTE(v) (uint8_t(Clamp(255.f * powf((v), 1.f/2.3f), 0.f, 255.f)))
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
}

/**
 * Loads a TGA image (not the RLE types though)
 */
 /*
static RGBSpectrum *ReadImageTGA(const string &name, int *width, int *height)
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
