/*
 * sampler.cpp
 * dnelson
 */
 
#include "sampler.h"
#include "integrator.h"
#include "volume.h"

// Sampler
Sampler::~Sampler() {
}

Sampler::Sampler(int xstart, int xend, int ystart, int yend, int spp, float sopen, float sclose)
    : xPixelStart(xstart), xPixelEnd(xend), yPixelStart(ystart),
      yPixelEnd(yend), samplesPerPixel(spp), shutterOpen(sopen),
      shutterClose(sclose)
{
		IF_DEBUG(cout << "Sampler(" << xstart << ", " << xend << ", " << ystart << ", " << yend << ", " << spp
									<< ", " << sopen << ", " << sclose << ") constructor." << endl);
}

bool Sampler::ReportResults(Sample *samples, const Ray *rays, const Spectrum *Ls, int count)
{
    return true;
}

void Sampler::ComputeSubWindow(int num, int count, int *newXStart,
        int *newXEnd, int *newYStart, int *newYEnd) const
{
    // Determine how many tiles to use in each dimension, _nx_ and _ny_
    int dx = xPixelEnd - xPixelStart, dy = yPixelEnd - yPixelStart;
    int nx = count, ny = 1;
		
    while ((nx & 0x1) == 0 && 2 * dx * ny < dy * nx) {
        nx >>= 1;
        ny <<= 1;
    }

    // Compute $x$ and $y$ pixel sample range for sub-window
    int xo = num % nx, yo = num / nx;
    float tx0 = float(xo) / float(nx), tx1 = float(xo+1) / float(nx);
    float ty0 = float(yo) / float(ny), ty1 = float(yo+1) / float(ny);
		
		IF_DEBUG(cout << "ComputeSubWindow(" << num << "," << count << ") nx = " << nx << " ny = "
									<< ny << " xo = " << xo << " yo = " << yo << " tx0 = " << tx0 << " tx1 = "
									<< tx1 << " ty0 = " << ty0 << " ty1 = " << ty1 << endl);
				 
    *newXStart = (int)floorf(Lerp(tx0, xPixelStart, xPixelEnd));
    *newXEnd   = (int)floorf(Lerp(tx1, xPixelStart, xPixelEnd));
    *newYStart = (int)floorf(Lerp(ty0, yPixelStart, yPixelEnd));
    *newYEnd   = (int)floorf(Lerp(ty1, yPixelStart, yPixelEnd));
}

// Sample 
Sample::Sample(Sampler *sampler, VolumeIntegrator *vol, const Scene *scene)
{
		IF_DEBUG(cout << "Sample() constructor." << endl);

    if (vol)  vol->RequestSamples(sampler, this, scene);
    AllocateSampleMemory();
}

void Sample::AllocateSampleMemory()
{
    // Allocate storage for sample pointers
    int nPtrs = n1D.size() + n2D.size();
    if (!nPtrs) {
        oneD = twoD = NULL;
        return;
    }
    oneD = AllocAligned<float *>(nPtrs);
    twoD = oneD + n1D.size();

    // Compute total number of sample values needed
    int totSamples = 0;
    for (uint32_t i = 0; i < n1D.size(); ++i)
        totSamples += n1D[i];
    for (uint32_t i = 0; i < n2D.size(); ++i)
        totSamples += 2 * n2D[i];
				
		IF_DEBUG(cout << "Sample::AllocateSampleMemory: nPtrs = " << nPtrs << " totSamples = " << totSamples << endl);

    // Allocate storage for sample values
    float *mem = AllocAligned<float>(totSamples);
    for (uint32_t i = 0; i < n1D.size(); ++i) {
        oneD[i] = mem;
        mem += n1D[i];
    }
    for (uint32_t i = 0; i < n2D.size(); ++i) {
        twoD[i] = mem;
        mem += 2 * n2D[i];
    }
}

Sample *Sample::Duplicate(int count) const
{
    Sample *ret = new Sample[count];
    for (int i = 0; i < count; ++i) {
        ret[i].n1D = n1D;
        ret[i].n2D = n2D;
        ret[i].AllocateSampleMemory();
    }
    return ret;
}

// StratifiedSampler
StratifiedSampler::StratifiedSampler(int xstart, int xend, int ystart, int yend, int xs, int ys, 
																		bool jitter, float sopen, float sclose)
    : Sampler(xstart, xend, ystart, yend, xs * ys, sopen, sclose)
{
		IF_DEBUG(cout << "StratifiedSampler(" << xstart << ", " << xend << ", " << ystart << ", " << yend << ", " << xs
									<< ", " << ys << ", " << jitter << ", " << sopen << ", " << sclose << ") constructor." << endl);

    jitterSamples = jitter;
    xPos = xPixelStart;
    yPos = yPixelStart;
    xPixelSamples = xs;
    yPixelSamples = ys;
    sampleBuf = new float[5 * xPixelSamples * yPixelSamples];
}


StratifiedSampler::~StratifiedSampler() {
    delete[] sampleBuf;
}


Sampler *StratifiedSampler::GetSubSampler(int num, int count)
{
    int x0, x1, y0, y1;
    ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
		
		IF_DEBUG(cout << "StratifiedSampler:GetSubSampler(" << num << "," << count << ") x0 = " << x0 << " x1 = " << x1
									<< " y0 = " << y0 << " y1 = " << y1 << endl);
		
    if (x0 == x1 || y0 == y1) return NULL;
		
    return new StratifiedSampler(x0, x1, y0, y1, xPixelSamples,
        yPixelSamples, jitterSamples, shutterOpen, shutterClose);
}


int StratifiedSampler::GetMoreSamples(Sample *samples, RNG &rng)
{
    if (yPos == yPixelEnd) return 0;
    int nSamples = xPixelSamples * yPixelSamples;
    // Generate stratified camera samples for _(xPos, yPos)_

    // Generate initial stratified samples into _sampleBuf_ memory
    float *bufp = sampleBuf;
    float *imageSamples = bufp; bufp += 2 * nSamples;
    float *lensSamples = bufp;  bufp += 2 * nSamples;
    float *timeSamples = bufp;
    StratifiedSample2D(imageSamples, xPixelSamples, yPixelSamples, rng, jitterSamples);
    StratifiedSample2D(lensSamples, xPixelSamples, yPixelSamples, rng, jitterSamples);
    StratifiedSample1D(timeSamples, xPixelSamples * yPixelSamples, rng, jitterSamples);

    // Shift stratified image samples to pixel coordinates
    for (int o = 0; o < 2 * xPixelSamples * yPixelSamples; o += 2) {
        imageSamples[o]   += xPos;
        imageSamples[o+1] += yPos;
    }

    // Decorrelate sample dimensions
    Shuffle(lensSamples, xPixelSamples*yPixelSamples, 2, rng);
    Shuffle(timeSamples, xPixelSamples*yPixelSamples, 1, rng);

    // Initialize stratified _samples_ with sample values
    for (int i = 0; i < nSamples; ++i) {
        samples[i].imageX = imageSamples[2*i];
        samples[i].imageY = imageSamples[2*i+1];
        samples[i].lensU = lensSamples[2*i];
        samples[i].lensV = lensSamples[2*i+1];
        samples[i].time = Lerp(timeSamples[i], shutterOpen, shutterClose);
				
				IF_DEBUG(cout << "GetMoreSamples (i=" << i << ") imageX=" << samples[i].imageX << " imageY=" 
											<< samples[i].imageY << " lensU=" << samples[i].lensU << " lensV=" << samples[i].lensV 
											<< " time=" << samples[i].time << endl);
				
        // Generate stratified samples for integrators
        for (uint32_t j = 0; j < samples[i].n1D.size(); ++j)
            LatinHypercube(samples[i].oneD[j], samples[i].n1D[j], 1, rng);
        for (uint32_t j = 0; j < samples[i].n2D.size(); ++j)
            LatinHypercube(samples[i].twoD[j], samples[i].n2D[j], 2, rng);
    }

    // Advance to next pixel for stratified sampling
    if (++xPos == xPixelEnd) {
        xPos = xPixelStart;
        ++yPos;
    }
    return nSamples;
}


StratifiedSampler *CreateStratifiedSampler(const Film *film, const Camera *camera)
{
    bool jitter = false;
		
    // Initialize common sampler parameters
    int xstart, xend, ystart, yend;
    film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
		
    int xsamp = 1; //2 samples per pixel
    int ysamp = 1; //2
		
    return new StratifiedSampler(xstart, xend, ystart, yend, xsamp, ysamp,
        jitter, camera->shutterOpen, camera->shutterClose);
}
