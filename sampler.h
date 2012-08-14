/*
 * sampler.h
 * dnelson
 */
 
#ifndef AREPO_RT_SAMPLER_H
#define AREPO_RT_SAMPLER_H

#include "ArepoRT.h"
#include "util.h"

class Sampler {
public:
    // construction
    virtual ~Sampler();
    Sampler(int xstart, int xend, int ystart, int yend, int spp, float sopen, float sclose);
    virtual int GetMoreSamples(Sample *sample, RNG &rng) = 0;
    virtual int MaximumSampleCount() = 0;
    virtual bool ReportResults(Sample *samples, const Ray *rays, const Spectrum *Ls, int count);
    virtual Sampler *GetSubSampler(int num, int count) = 0;
    virtual int RoundSize(int size) const = 0;

    // data
    const int xPixelStart, xPixelEnd, yPixelStart, yPixelEnd;
    const int samplesPerPixel;
    const float shutterOpen, shutterClose;
protected:
    void ComputeSubWindow(int num, int count, int *xstart, int *xend, int *ystart, int *yend) const;
};

struct CameraSample {
    float imageX, imageY;
    float lensU, lensV;
    float time;
};

struct Sample : public CameraSample {
    // public methods
    Sample(Sampler *sampler, VolumeIntegrator *vol, const Scene *scene);
    uint32_t Add1D(uint32_t num) {
        n1D.push_back(num);
        return n1D.size()-1;
    }
    uint32_t Add2D(uint32_t num) {
        n2D.push_back(num);
        return n2D.size()-1;
    }
    ~Sample() {
        if (oneD != NULL) {
            FreeAligned(oneD[0]);
            FreeAligned(oneD);
        }
    }
    Sample *Duplicate(int count) const;

    // data
    vector<uint32_t> n1D, n2D;
    float **oneD, **twoD;
		
private:
    // private methods
    void AllocateSampleMemory();
    Sample() { oneD = twoD = NULL; }
};

class StratifiedSampler : public Sampler {
public:
    // construction
    StratifiedSampler(int xstart, int xend, int ystart, int yend,
                      int xs, int ys, bool jitter, float sopen, float sclose);
    ~StratifiedSampler();
		
		// methods
    int RoundSize(int size) const { return size; }
    Sampler *GetSubSampler(int num, int count);
    int GetMoreSamples(Sample *sample, RNG &rng);
    int MaximumSampleCount() { return xPixelSamples * yPixelSamples; }
private:
    // data
    int xPixelSamples, yPixelSamples;
    bool jitterSamples;
    int xPos, yPos;
    float *sampleBuf;
};

StratifiedSampler *CreateStratifiedSampler(const Film *film, const Camera *camera);

#endif
