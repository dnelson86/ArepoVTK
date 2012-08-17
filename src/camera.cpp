/*
 * camera.cpp
 * dnelson
 */

#include "ArepoRT.h"
#include "geometry.h"
#include "transform.h"
#include "sampler.h"
#include "fileio.h" 
#include "camera.h"
#include "spectrum.h"

// Filter

Filter::~Filter() {
}

float BoxFilter::Evaluate(float x, float y) const
{
    return 1.0;
}

BoxFilter *CreateBoxFilter()
{
    float xw = 0.5f;
    float yw = 0.5f;
    return new BoxFilter(xw, yw);
}

// Film

Film::Film(int xres, int yres, Filter *filt, const float crop[4], const string &fn, bool openWindow)
    : xResolution(xres), yResolution(yres)
{
		IF_DEBUG(cout << "Film(" << xres << ", " << yres << ", ...) constructor." << endl);
		
    filter = filt;
    memcpy(cropWindow, crop, 4 * sizeof(float));
    filename = fn;
		
    // Compute film image extent
    xPixelStart = (int)ceilf(xResolution * cropWindow[0]);
    xPixelCount = max(1, (int)ceilf(xResolution * cropWindow[1]) - xPixelStart);
    yPixelStart = (int)ceilf(yResolution * cropWindow[2]);
    yPixelCount = max(1, (int)ceilf(yResolution * cropWindow[3]) - yPixelStart);

		IF_DEBUG(cout << " xPixelStart = " << xPixelStart << " xPixelCount = " << xPixelCount
									<< " yPixelStart = " << yPixelStart << " yPixelCount = " << yPixelCount << endl);
		
    // Allocate film image storage
    pixels = new BlockedArray<Pixel>(xPixelCount, yPixelCount);

    // Precompute filter weight table
    filterTable = new float[FILTER_TABLE_SIZE * FILTER_TABLE_SIZE];
    float *ftp = filterTable;
		
    for (int y = 0; y < FILTER_TABLE_SIZE; ++y)
		{
        float fy = ((float)y + 0.5f) * filter->yWidth / FILTER_TABLE_SIZE;
        for (int x = 0; x < FILTER_TABLE_SIZE; ++x)
				{
            float fx = ((float)x + 0.5f) * filter->xWidth / FILTER_TABLE_SIZE;
            *ftp++ = filter->Evaluate(fx, fy);
        }
    }

    // Possibly open window for image display
    if (openWindow)
		{
        cout << "WARNING: Window support not yet implemented." << endl;
    }
}


void Film::AddSample(const CameraSample &sample, const Spectrum &L)
{				 
    // Compute sample's raster extent
    float dimageX = sample.imageX - 0.5f;
    float dimageY = sample.imageY - 0.5f;
		
    int x0 = (int)ceilf(dimageX - filter->xWidth);
    int x1 = (int)floorf(dimageX + filter->xWidth);
    int y0 = (int)ceilf(dimageY - filter->yWidth);
    int y1 = (int)floorf(dimageY + filter->yWidth);
		
		IF_DEBUG(cout << "Film::AddSample(imageX=" << sample.imageX << " imageY=" << sample.imageY 
									<< ") (x0=" << x0 << " x1=" << x1 << " y0=" << y0
									<< " y1=" << y1 << ")" << endl);
		
    x0 = max(x0, xPixelStart);
    x1 = min(x1, xPixelStart + xPixelCount - 1);
    y0 = max(y0, yPixelStart);
    y1 = min(y1, yPixelStart + yPixelCount - 1);
		
		IF_DEBUG(cout << " (xPixelStart=" << xPixelStart << " xPixelCount=" << xPixelCount
									<< " (yPixelStart=" << yPixelStart << " yPixelCount=" << yPixelCount
									<< ") (x0=" << x0 << " x1=" << x1 << " y0=" << y0
									<< " y1=" << y1 << ")" << endl);
		
    if ((x1-x0) < 0 || (y1-y0) < 0)
    {
				IF_DEBUG(cout << " WARNING: Sample outside image extent." << endl);
        return;
    }

    // Loop over filter support and add sample to pixel arrays
    float xyz[3];
    L.ToXYZ(xyz);

    // Precompute $x$ and $y$ filter table offsets
    int *ifx = (int *)alloca((x1 - x0 + 1) * sizeof(int));
    for (int x = x0; x <= x1; ++x)
		{
        float fx = fabsf((x - dimageX) * filter->invXWidth * FILTER_TABLE_SIZE);
        ifx[x-x0] = min((int)floorf(fx), FILTER_TABLE_SIZE-1);
    }
		
    int *ify = (int *)alloca((y1 - y0 + 1) * sizeof(int));
    for (int y = y0; y <= y1; ++y)
		{
        float fy = fabsf((y - dimageY) * filter->invYWidth * FILTER_TABLE_SIZE);
        ify[y-y0] = min((int)floorf(fy), FILTER_TABLE_SIZE-1);
    }
		
    bool syncNeeded = (filter->xWidth > 0.5f || filter->yWidth > 0.5f);
		
    for (int y = y0; y <= y1; ++y)
		{
        for (int x = x0; x <= x1; ++x)
				{
            // Evaluate filter value at $(x,y)$ pixel
            int offset = ify[y-y0] * FILTER_TABLE_SIZE + ifx[x-x0];
            float filterWt = filterTable[offset];

            // Update pixel values with filtered sample contribution
            Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);
            if (!syncNeeded) {
                pixel.Lxyz[0] += filterWt * xyz[0];
                pixel.Lxyz[1] += filterWt * xyz[1];
                pixel.Lxyz[2] += filterWt * xyz[2];
                pixel.weightSum += filterWt;
            }
            else {
                // Safely update _Lxyz_ and _weightSum_ even with concurrency
                //AtomicAdd(&pixel.Lxyz[0], filterWt * xyz[0]);
                //AtomicAdd(&pixel.Lxyz[1], filterWt * xyz[1]);
                //AtomicAdd(&pixel.Lxyz[2], filterWt * xyz[2]);
                //AtomicAdd(&pixel.weightSum, filterWt);
								cout << "WARNING: AtomicAdd1 not implemented." << endl;
            }
        }
    }
}

bool Film::DrawLine(float x1, float y1, float x2, float y2, const Spectrum &L)
{
		// (x1,y1) and (x2,y2) line segment endpoints in raster space
		IF_DEBUG(cout << "Film::DrawLine(" << x1 << "," << y1 << "," << x2 << "," << y2 << ")" << endl);
		
		// convert spectrum to XYZ and set flat filter
    float xyz[3];
    L.ToXYZ(xyz);
		
		float filterWt = 1.0;
		
		// premultiplies
		xyz[0] *= filterWt;
		xyz[1] *= filterWt;
		xyz[2] *= filterWt;
		
		//bool syncNeeded = false; // not thread safe!
		int off1, off2;	

#ifdef USE_LINEALGO_BRESENHAM
		// algorithm: Bresenham, Jack E. "Algorithm for computer control of a digital plotter".
		//					  IBM Systems Journal, Vol. 4, No.1, January 1965, pp 25-30.	
		
		const bool steep = (fabs(y2-y1) > fabs(x2-x1));
		if (steep) {
				swap(x1,y1);
				swap(x2,y2);
		}
		if (x1 > x2) {
				swap(x1,x2);
				swap(y1,y2);
		}
		
		const float dx = x2 - x1;
		const float dy = fabs(y2-y1);
		
		float error = dx / 2.0f;
		const int ystep = (y1 < y2) ? 1 : -1;
		int y = (int)y1;
		
		const int max_x = (int)x2;
		const int max_y = (int)y2;
	
		for (int x=(int)x1; x < max_x; x++) {
				IF_DEBUG(cout << " Painting x = " << x << " y = " << y << endl);
				
				if(steep) {
						off1 = (y - yPixelStart);	//(y,x)
						off2 = (x - xPixelStart);
				}
				else {
						off1 = (x - xPixelStart);	//(x,y)
						off2 = (y - yPixelStart);
				}
				
				if (off1 < 0 || off2 < 0 || off1 >= xPixelCount || off2 >= yPixelCount) {
						IF_DEBUG(cout << "WARNING: Painting line outside image extent (" 
						              << off1 << "," << off2 << ")" << endl);
				} else {
						Pixel &pixel = (*pixels)(off1,off2);

						if (!syncNeeded) {
								IF_DEBUG(cout << " pixel Lxyz before: " << pixel.Lxyz[0] << " " << pixel.Lxyz[1] << " " 
															<< pixel.Lxyz[2] << endl);
								pixel.Lxyz[0] += xyz[0];
								pixel.Lxyz[1] += xyz[1];
								pixel.Lxyz[2] += xyz[2];
								IF_DEBUG(cout << " pixel Lxyz after:  " << pixel.Lxyz[0] << " " << pixel.Lxyz[1] << " " 
															<< pixel.Lxyz[2] << endl);
								pixel.weightSum += filterWt;
						}
				}
				
				error -= dy;
				if (error < 0) {
						y += ystep;
						error += dx;
				}
		}
#endif
#ifdef USE_LINEALGO_WU

		float dx = x2 - x1;
		float dy = y2 - y1;
		const bool steep = (fabs(dy) > fabs(dx));
		
		// special case for line extent of 1 pixel or less
		if ( fabs(dx) < 1.0f && fabs(dy) < 1.0f ) {
				IF_DEBUG(cout << " WARNING: Line extent sub-pixel." << endl);
        return false;
		}

		// algorithm: Wu, Xiaolin (July 1991): "An efficient antialiasing technique".
		//            Computer Graphics 25 (4): 143-152. doi:10.1145/127719.122734.
		// note: should handle vertical, horizontal, and 45 deg lines as special cases (no need for AA)
		
		if ( steep ) {
				swap(x1,y1);
				swap(x2,y2);
				swap(dx,dy);
		}
		
		if ( x2 < x1 ) {
				swap(x1,x2);
				swap(y1,y2);
		}
		
#define ipart_(x)  ((int)(x))
#define round_(x)  ((int)(((float)(x))+0.5f))
#define fpart_(x)  (((float)(x))-(float)ipart_(x))
#define rfpart_(x) (1.0f-fpart_(x))		
		
		float gradient = dy / dx;		
		float xend,yend,xgap;//,ygap;
		
/*	// special case: vertical line (dy==0 because swapped already)
		IF_DEBUG(cout << " dx: " << dx << " dy: " << dy << " grad: " << gradient << endl);
		if (dy == 0 || fabs(dy) < 1e-4) {
				IF_DEBUG(cout << " VERTICAL " << (int)x1 << " " << (int)x2 << " " << round_(y1) << endl);
				for (int x=(int)x1; x<(int)x2; x++) {
						if (round_(y1) > 0 && x > 0 && round_(y1) < xPixelCount && x < yPixelCount) {
						    IF_DEBUG(cout << "  Paint " << x << " " << round_(y1) << endl);
								Pixel &pixel = (*pixels)(round_(y1),x);
								pixel.Lxyz[0] += filterWt * xyz[0];
								pixel.Lxyz[1] += filterWt * xyz[1] * 1.0;
								pixel.Lxyz[2] += filterWt * xyz[2];
								pixel.weightSum += filterWt;
						}
				}
				return true;
		} */
		
		// first endpoint
		xend = round_(x1);
		yend = y1 + gradient * (xend-x1);
		xgap = rfpart_(x1 + 0.5f);

		int xpx11 = (int)xend;
		int ypx11 = ipart_(yend);
		
		off1 = xpx11;	//(y,x)
		off2 = ypx11;

		if (steep) swap(off1,off2); //(x,y)		
		
		if (off1 > 0 && off2 > 0 && off1 < xPixelCount && off2 < yPixelCount) {		
				Pixel &pixel = (*pixels)(off1,off2);
				pixel.Lxyz[0] += xyz[0] * rfpart_(yend) * xgap;
				// multiplicative brightness modification on y-comp as per CIE XYZ standard
				// note: apparently not, changes colors (not strict brightness modification)
				//       better to multiply all components
				pixel.Lxyz[1] += xyz[1] * rfpart_(yend) * xgap;
				pixel.Lxyz[2] += xyz[2] * rfpart_(yend) * xgap;
				pixel.weightSum += filterWt;
		}
		if (off1 > 0 && off2+1 > 0 && off1 < xPixelCount && off2+1 < yPixelCount) {		
				Pixel &pixel = (*pixels)(off1,off2+1);
				pixel.Lxyz[0] += xyz[0] * fpart_(yend) * xgap;
				pixel.Lxyz[1] += xyz[1] * fpart_(yend) * xgap;
				pixel.Lxyz[2] += xyz[2] * fpart_(yend) * xgap;
				pixel.weightSum += filterWt;
		}
		
		float intery = yend + gradient;
		
		// second endpoint
		xend = round_(x2);
		yend = y2 + gradient * (xend-x2);
		xgap = fpart_(x2 + 0.5f);
		
		int xpx12 = (int)xend;
		int ypx12 = ipart_(yend);
		
		off1 = xpx12;	//(y,x)
		off2 = ypx12;

		if (steep) swap(off1,off2); //(x,y)	
		
		if (off1 > 0 && off2 > 0 && off1 < xPixelCount && off2 < yPixelCount) {		
				Pixel &pixel = (*pixels)(off1,off2);
				pixel.Lxyz[0] += xyz[0] * rfpart_(yend) * xgap;
				pixel.Lxyz[1] += xyz[1] * rfpart_(yend) * xgap;
				pixel.Lxyz[2] += xyz[2] * rfpart_(yend) * xgap;
				pixel.weightSum += filterWt;
		}
		if (off1 > 0 && off2+1 > 0 && off1 < xPixelCount && off2+1 < yPixelCount) {		
				Pixel &pixel = (*pixels)(off1,off2+1);
				pixel.Lxyz[0] += xyz[0] * fpart_(yend) * xgap;
				pixel.Lxyz[1] += xyz[1] * fpart_(yend) * xgap;
				pixel.Lxyz[2] += xyz[2] * fpart_(yend) * xgap;
				pixel.weightSum += filterWt;
		}
		
		// main loop
		for (int x = xpx11+1; x <= (xpx12-1); x++)
		{
				off1 = x;	//(y,x)
				off2 = ipart_(intery);

				if (steep) swap(off1,off2); //(x,y)

				if (off1 > 0 && off2 > 0 && off1 < xPixelCount && off2 < yPixelCount) {		
						Pixel &pixel = (*pixels)(off1,off2);
						pixel.Lxyz[0] += xyz[0] * rfpart_(intery);
						pixel.Lxyz[1] += xyz[1] * rfpart_(intery);
						pixel.Lxyz[2] += xyz[2] * rfpart_(intery);
						pixel.weightSum = filterWt;
						IF_DEBUG(cout << " Painting x = " << off1 << " y = " << off2 << " (" 
						              << rfpart_(intery) << ")" << endl);
				}
				if (off1 > 0 && off2+1 > 0 && off1 < xPixelCount && off2+1 < yPixelCount) {	
						Pixel &pixel = (*pixels)(off1,off2+1);
						pixel.Lxyz[0] += xyz[0] * fpart_(intery);
						pixel.Lxyz[1] += xyz[1] * fpart_(intery);
						pixel.Lxyz[2] += xyz[2] * fpart_(intery);
						pixel.weightSum = filterWt;
						IF_DEBUG(cout << " Painting x = " << off1 << " y = " << off2+1 << " (" 
						              << fpart_(intery) << ")" << endl);
				}
				intery += gradient;
		}
		
#endif
		
		return true;
}

void Film::Splat(const CameraSample &sample, const Spectrum &L)
{
		// NOT IMPLEMENTED
    float xyz[3];
    L.ToXYZ(xyz);
		
    int x = (int)floorf(sample.imageX);
		int y = (int)floorf(sample.imageY);
    if (x < xPixelStart || x - xPixelStart >= xPixelCount ||
        y < yPixelStart || y - yPixelStart >= yPixelCount) return;
				
    //Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);
		
    //AtomicAdd(&pixel.splatXYZ[0], xyz[0]);
    //AtomicAdd(&pixel.splatXYZ[1], xyz[1]);
    //AtomicAdd(&pixel.splatXYZ[2], xyz[2]);
		cout << "WARNING: AtomicAdd2 not yet." << endl;
}

void Film::GetSampleExtent(int *xstart, int *xend, int *ystart, int *yend) const
{
    *xstart = (int)floorf(xPixelStart + 0.5f - filter->xWidth);
    //*xend   = (int)floorf(xPixelStart + 0.5f + xPixelCount  + filter->xWidth); //bug? too big by 1 (dnelson)
		*xend   = (int)floorf(xPixelStart + 0.5f + xPixelCount-1  + filter->xWidth);

    *ystart = (int)floorf(yPixelStart + 0.5f - filter->yWidth);
    //*yend   = (int)floorf(yPixelStart + 0.5f + yPixelCount + filter->yWidth); //bug? too big by 1 (dnelson)
		*yend   = (int)floorf(yPixelStart + 0.5f + yPixelCount-1 + filter->yWidth);
}

void Film::GetPixelExtent(int *xstart, int *xend, int *ystart, int *yend) const
{
    *xstart = xPixelStart;
    *xend   = xPixelStart + xPixelCount;
    *ystart = yPixelStart;
    *yend   = yPixelStart + yPixelCount;
}


void Film::WriteImage(float splatScale)
{
		IF_DEBUG(cout << "Film:WriteImage(" << splatScale << ") nx = " << xPixelCount << " ny = " 
		              << yPixelCount << endl);
		
    // Convert image to RGB and compute final pixel values
    int nPix = xPixelCount * yPixelCount;
    float *rgb = new float[3*nPix];
    int offset = 0;
		
    for (int y = 0; y < yPixelCount; ++y) {
        for (int x = 0; x < xPixelCount; ++x) {
            // Convert pixel XYZ color to RGB
            XYZToRGB((*pixels)(x, y).Lxyz, &rgb[3*offset]);

            // Normalize pixel with weight sum
            float weightSum = (*pixels)(x, y).weightSum;
            if (weightSum != 0.0f) {
                float invWt = 1.0f / weightSum;
                rgb[3*offset  ] = max(0.0f, rgb[3*offset  ] * invWt);
                rgb[3*offset+1] = max(0.0f, rgb[3*offset+1] * invWt);
                rgb[3*offset+2] = max(0.0f, rgb[3*offset+2] * invWt);
            }

            // Add splat value at pixel
            float splatRGB[3];
            XYZToRGB((*pixels)(x, y).splatXYZ, splatRGB);
            rgb[3*offset  ] += splatScale * splatRGB[0];
            rgb[3*offset+1] += splatScale * splatRGB[1];
            rgb[3*offset+2] += splatScale * splatRGB[2];
            ++offset;
        }
    }

    // Write RGB image
    ::WriteImage(filename, rgb, NULL, xPixelCount, yPixelCount,
                 xResolution, yResolution, xPixelStart, yPixelStart);

    // Release temporary image memory
    delete[] rgb;
}

void Film::WriteRawRGB()
{
		IF_DEBUG(cout << "Film:WriteRawRGB()" << endl);
		
		// Convert image to RGB and compute final pixel values
    float *rgb = new float[3*xPixelCount*yPixelCount];
		int offset = 0;
		
    for (int y = 0; y < yPixelCount; ++y) {
        for (int x = 0; x < xPixelCount; ++x) {
            // Convert pixel XYZ color to RGB
            XYZToRGB((*pixels)(x, y).Lxyz, &rgb[3*offset]);

            // Normalize pixel with weight sum
            float weightSum = (*pixels)(x, y).weightSum;
            if (weightSum != 0.0f) {
                float invWt = 1.0f / weightSum;
                rgb[3*offset  ] = max(0.0f, rgb[3*offset  ] * invWt);
                rgb[3*offset+1] = max(0.0f, rgb[3*offset+1] * invWt);
                rgb[3*offset+2] = max(0.0f, rgb[3*offset+2] * invWt);
            }

            ++offset;
        }
    }
		
		// Write RAW text file
		WriteFloatFile(Config.rawRGBFile.c_str(), rgb, xPixelCount, yPixelCount);
		
		// Release temporary image memory
		delete[] rgb;

}

void Film::UpdateDisplay(int x0, int y0, int x1, int y1, float splatScale)
{
		IF_DEBUG(cout << "Film::UpdateDisplay()" << endl);
}

Film *CreateFilm(Filter *filter)
{
    string filename = Config.imageFile;
    if (filename == "")
#ifdef HAVE_OPENEXR
        filename = "frame.exr";
#else
        filename = "frame.tga";
#endif

    int xres = Config.imageXPixels;
    int yres = Config.imageYPixels;
    bool openwin = Config.openWindow;
		float crop[4] = { 0, 1, 0, 1 }; // get from Config
		
    return new Film(xres, yres, filter, crop, filename, openwin);
}

// Camera
Camera::~Camera() {
    delete film;
}

Camera::Camera(const Transform &cam2world, float sopen, float sclose, Film *f)
    : CameraToWorld(cam2world), shutterOpen(sopen), shutterClose(sclose)
{
		IF_DEBUG(cout << "Camera(c2w," << sopen << "," << sclose << ",f) constructor." << endl);
		
    film = f;
    //if (CameraToWorld.HasScale())
    //    cout << "Warning! CameraToWorld has a scale factor." << endl;
}

bool Camera::RasterizeLine(const Point &p1, const Point &p2, const Spectrum &L)
{
		p1.print("Camera::RasterizeLine p1 W ");
		p2.print("Camera::RasterizeLine p2 W ");
		
		// transform start and end points to raster space
    //float x1,y1,x2,y2;
		Point start,end;
		
		Transform w2r = Inverse(RasterToCamera) * Inverse(CameraToWorld);
		start = w2r(p1);
		end   = w2r(p2);
		
		w2r.print(" w2r ");
		Inverse(CameraToWorld).print(" inv c2w ");
		Inverse(RasterToCamera).print(" inc r2c ");
		
		start.print("Camera::RasterizeLine p1 R ");
		end.print("Camera::RasterizeLine p2 R ");

		// ask film to draw
		return film->DrawLine(start.x,start.y,end.x,end.y,L);
}

Camera::Camera(const Transform &cam2world, const Transform &proj, const float screenWindow[4],
               float sopen, float sclose, float lensr, float focald, Film *f)
		: CameraToWorld(cam2world), shutterOpen(sopen), shutterClose(sclose)
{
		IF_DEBUG(cout << "Camera(c2w,proj,sW," << sopen << "," << sclose << "," << lensr << "," << focald 
		              << "f) constructor." << endl);
		
		film = f;
	
    // Initialize depth of field parameters
    lensRadius = lensr;
    focalDistance = focald;

    // Compute projective camera transformations
    CameraToScreen = proj;

    // Compute projective camera screen transformations
    ScreenToRaster = Scale(float(film->xResolution),
                           float(film->yResolution), 1.0f) *
        Scale(1.0f / (screenWindow[1] - screenWindow[0]),
              1.0f / (screenWindow[2] - screenWindow[3]), 1.0f) *
        Translate(Vector(-screenWindow[0], -screenWindow[3], 0.0f));
		ScreenToRaster.print(" ScreenToRaster ");
    RasterToScreen = Inverse(ScreenToRaster);
    RasterToCamera = Inverse(CameraToScreen) * RasterToScreen;
}

ProjCamera::ProjCamera(const Transform &cam2world, const Transform &proj, const float screenWindow[4],
										   float sopen, float sclose, float lensr, float focald, Film *f)
		: Camera(cam2world, sopen, sclose, f)
{
		IF_DEBUG(cout << "ProjCamera(c2w,proj,sW," << sopen << "," << sclose << "," << lensr << "," << focald 
		              << "f) constructor." << endl);
		
    // Initialize depth of field parameters
    lensRadius = lensr;
    focalDistance = focald;

    // Compute projective camera transformations
    CameraToScreen = proj;

    // Compute projective camera screen transformations
    ScreenToRaster = Scale(float(film->xResolution),
                           float(film->yResolution), 1.0f) *
        Scale(1.0f / (screenWindow[1] - screenWindow[0]),
              1.0f / (screenWindow[2] - screenWindow[3]), 1.0f) *
        Translate(Vector(-screenWindow[0], -screenWindow[3], 0.0f));
		ScreenToRaster.print(" ScreenToRaster ");
    RasterToScreen = Inverse(ScreenToRaster);
    RasterToCamera = Inverse(CameraToScreen) * RasterToScreen;
}

// OrthoCamera
OrthoCamera::OrthoCamera(const Transform &cam2world, const float screenWindow[4], 
												 float sopen, float sclose, float lensr, float focald, Film *f)
    : Camera(cam2world, Orthographic(0.0, 1.0), screenWindow, sopen, sclose, lensr, focald, f)
{
		IF_DEBUG(cout << "OrthoCamera(c2w,sW=[" << screenWindow[0] << "," << screenWindow[1] << ","
	                << screenWindow[2] << "," << screenWindow[3] << "],...) constructor." << endl);

    // Compute differential changes in origin for ortho camera rays
    dxCamera = RasterToCamera(Vector(1, 0, 0));
    dyCamera = RasterToCamera(Vector(0, 1, 0));
}


float OrthoCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{
		IF_DEBUG(cout << "OrthoCamera::GenerateRay()" << endl);
		
    // Generate raster and camera samples
    Point Pras(sample.imageX, sample.imageY, 0);
    Point Pcamera;
		
    RasterToCamera(Pras, &Pcamera);
		IF_DEBUG(Pras.print(" sample (raster): "));
		IF_DEBUG(Pcamera.print(" sample (camera): "));
    *ray = Ray(Pcamera, Vector(0,0,1), 0.0f, INFINITY);
		
    if (lensRadius > 0.0) {
				// Modify ray for depth of field
				// NOT IMPLEMENTED
				cout << "ERROR: OrthoCamera GenerateRay lensRadius not zero." << endl;
    }
		
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
		
		IF_DEBUG(ray->printRay(" genray C "));
		
    CameraToWorld(*ray, ray);
				
		IF_DEBUG(ray->printRay(" genray W "));
				 
    return 1.0f;
}

// PerspectiveCamera
PerspectiveCamera::PerspectiveCamera(const Transform &cam2world, const float screenWindow[4], 
										float sopen, float sclose, float lensr, float focald, float fov, Film *film)
    : ProjCamera(cam2world, Perspective(fov, 1e-2f, 1000.0f), screenWindow,
								 sopen, sclose, lensr, focald, film)
{
		IF_DEBUG(cout << "PerspectiveCamera(c2w,sW=[" << screenWindow[0] << "," << screenWindow[1] << ","
	                << screenWindow[2] << "," << screenWindow[3] << "],...) constructor." << endl);

    // Compute differential changes in origin for ortho camera rays
    dxCamera = RasterToCamera(Point(1,0,0)) - RasterToCamera(Point(0,0,0));
    dyCamera = RasterToCamera(Point(0,1,0)) - RasterToCamera(Point(0,0,0));
}

float PerspectiveCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{
		IF_DEBUG(cout << "PerspectiveCamera::GenerateRay()" << endl);
		
    // Generate raster and camera samples
    Point Pras(sample.imageX, sample.imageY, 0);
    Point Pcamera;
    RasterToCamera(Pras, &Pcamera);
		IF_DEBUG(Pras.print(" sample (raster): "));
		IF_DEBUG(Pcamera.print(" sample (camera): "));
		
    *ray = Ray(Point(0,0,0), Normalize(Vector(Pcamera)), 0.0f, INFINITY);

    if (lensRadius > 0.0) {
				// Modify ray for depth of field		
				// NOT IMPLEMENTED
				cout << "ERROR: PerspectiveCamera GenerateRay lensRadius not zero." << endl;
    }
		
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
		
		IF_DEBUG(ray->printRay(" genray C "));
		
    CameraToWorld(*ray, ray);
				
		IF_DEBUG(ray->printRay(" genray W "));

    return 1.0f;
}

// Camera object creation
OrthoCamera *CreateOrthographicCamera(const Transform &cam2world, Film *film)
{
		// Config
    float shutteropen = 0.0f;
    float shutterclose = 1.0f;

    float lensradius = 0.0f;
    float focaldistance = 1e30f;
    float frame = float(film->xResolution)/float(film->yResolution);
    float screen[4];
		
    if (frame > 1.0f) {
        screen[0] = -frame;
        screen[1] =  frame;
        screen[2] = -1.0f;
        screen[3] =  1.0f;
    }
    else {
        screen[0] = -1.0f;
        screen[1] =  1.0f;
        screen[2] = -1.0f / frame;
        screen[3] =  1.0f / frame;
    }
		
    for (int i=0; i < 4; i++)
		  screen[i] *= Config.swScale;
				
    return new OrthoCamera(cam2world, screen, shutteropen, shutterclose,
        lensradius, focaldistance, film);
}

PerspectiveCamera *CreatePerspectiveCamera(const Transform &cam2world, Film *film)
{
    // Config
    float shutteropen = 0.0f;
    float shutterclose = 1.0f;

    float lensradius = 0.0f;
    float focaldistance = 1e30f;
    float frame = float(film->xResolution)/float(film->yResolution);
    float screen[4];
		
    if (frame > 1.0f) {
        screen[0] = -frame;
        screen[1] =  frame;
        screen[2] = -1.0f;
        screen[3] =  1.0f;
    }
    else {
        screen[0] = -1.0f;
        screen[1] =  1.0f;
        screen[2] = -1.0f / frame;
        screen[3] =  1.0f / frame;
    }
		
    for (int i=0; i < 4; i++)
		  screen[i] *= Config.swScale;
			
    float fov = Config.cameraFOV;
		
    return new PerspectiveCamera(cam2world, screen, shutteropen, shutterclose, lensradius, 
																 focaldistance, fov, film);
}
