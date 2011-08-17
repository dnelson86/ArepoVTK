/*
 * util.h
 * dnelson
 */
 
#ifndef AREPO_RT_UTIL_H
#define AREPO_RT_UTIL_H

#include "ArepoRT.h"

// timing

class Timer {
public:
    // construction
    Timer();
    
		// methods
    void Start();
    void Stop();
    void Reset();
    
    double Time();
private:
    // data
    double time0, elapsed;
    bool running;
    double GetTime();
    // UNIX Timer
    struct timeval timeofday;
};

// temporary memory management

void *AllocAligned(size_t size);

template <typename T> T *AllocAligned(uint32_t count) {
    return (T *)AllocAligned(count * sizeof(T));
}

void FreeAligned(void *);

template <typename T, int logBlockSize> class BlockedArray {
public:
    // BlockedArray Public Methods
    BlockedArray(uint32_t nu, uint32_t nv, const T *d = NULL) {
        uRes = nu;
        vRes = nv;
        uBlocks = RoundUp(uRes) >> logBlockSize;
        uint32_t nAlloc = RoundUp(uRes) * RoundUp(vRes);
        data = AllocAligned<T>(nAlloc);
        for (uint32_t i = 0; i < nAlloc; ++i)
            new (&data[i]) T();
        if (d)
            for (uint32_t v = 0; v < vRes; ++v)
                for (uint32_t u = 0; u < uRes; ++u)
                    (*this)(u, v) = d[v * uRes + u];
    }
    uint32_t BlockSize() const { return 1 << logBlockSize; }
    uint32_t RoundUp(uint32_t x) const {
        return (x + BlockSize() - 1) & ~(BlockSize() - 1);
    }
    uint32_t uSize() const { return uRes; }
    uint32_t vSize() const { return vRes; }
    ~BlockedArray() {
        for (uint32_t i = 0; i < uRes * vRes; ++i)
            data[i].~T();
        FreeAligned(data);
    }
    uint32_t Block(uint32_t a) const { return a >> logBlockSize; }
    uint32_t Offset(uint32_t a) const { return (a & (BlockSize() - 1)); }
    T &operator()(uint32_t u, uint32_t v) {
        uint32_t bu = Block(u), bv = Block(v);
        uint32_t ou = Offset(u), ov = Offset(v);
        uint32_t offset = BlockSize() * BlockSize() * (uBlocks * bv + bu);
        offset += BlockSize() * ov + ou;
        return data[offset];
    }
    const T &operator()(uint32_t u, uint32_t v) const {
        uint32_t bu = Block(u), bv = Block(v);
        uint32_t ou = Offset(u), ov = Offset(v);
        uint32_t offset = BlockSize() * BlockSize() * (uBlocks * bv + bu);
        offset += BlockSize() * ov + ou;
        return data[offset];
    }
    void GetLinearArray(T *a) const {
        for (uint32_t v = 0; v < vRes; ++v)
            for (uint32_t u = 0; u < uRes; ++u)
                *a++ = (*this)(u, v);
    }
private:
    // BlockedArray Private Data
    T *data;
    uint32_t uRes, vRes, uBlocks;
};

// random number generation

class RNG {
public:
    RNG(uint32_t seed = 5489UL) {
				IF_DEBUG(cout << "RNG(" << seed << ") constructor." << endl);
        mti = N+1; /* mti==N+1 means mt[N] is not initialized */
        Seed(seed);
    }

    void Seed(uint32_t seed) const;
    float RandomFloat() const;
    unsigned long RandomUInt() const;

private:
    static const int N = 624;
    mutable unsigned long mt[N]; /* the array for the state vector  */
    mutable int mti;
};

static const float OneMinusEpsilon=0x1.fffffep-1;

void StratifiedSample1D(float *samples, int nsamples, RNG &rng, bool jitter = true);
void StratifiedSample2D(float *samples, int nx, int ny, RNG &rng, bool jitter = true);
void LatinHypercube(float *samples, uint32_t nSamples, uint32_t nDim, RNG &rng);

template <typename T> void Shuffle(T *samp, uint32_t count, uint32_t dims, RNG &rng)
{
    for (uint32_t i = 0; i < count; ++i) {
        uint32_t other = i + (rng.RandomUInt() % (count - i));
        for (uint32_t j = 0; j < dims; ++j)
            swap(samp[dims*i + j], samp[dims*other + j]);
    }
}


#endif
