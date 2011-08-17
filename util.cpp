/*
 * util.cpp
 * dnelson
 */
 
#include "util.h"

// Timer
Timer::Timer()
{
    time0 = elapsed = 0;
    running = 0;
}

double Timer::GetTime()
{
    gettimeofday( &timeofday, NULL );
    return timeofday.tv_sec + timeofday.tv_usec / 1000000.0;
}

void Timer::Start()
{
    running = 1;
    time0 = GetTime();
}

void Timer::Stop()
{
    running = 0;
    elapsed += GetTime() - time0;
}

void Timer::Reset()
{
    running = 0;
    elapsed = 0;
}

double Timer::Time()
{
    if (running) {
        Stop();
        Start();
    }
    return elapsed;
}

// memory

void *AllocAligned(size_t size) {
    return memalign(L1_CACHE_LINE_SIZE, size);
}

void FreeAligned(void *ptr) {
    if (!ptr) return;
    free(ptr);
}

// random number generation

#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

void RNG::Seed(uint32_t seed) const {
    mt[0]= seed & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
        (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* generates a random number on [0,1)-real-interval */
float RNG::RandomFloat() const
{
    float v = (RandomUInt() & 0xffffff) / float(1 << 24);
    return v;
}

unsigned long RNG::RandomUInt() const
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if Seed() has not been called, */
            Seed(5489UL); /* default initial seed */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

void StratifiedSample1D(float *samp, int nSamples, RNG &rng, bool jitter) {
    float invTot = 1.f / nSamples;
    for (int i = 0;  i < nSamples; ++i) {
        float delta = jitter ? rng.RandomFloat() : 0.5f;
        *samp++ = min((i + delta) * invTot, OneMinusEpsilon);
    }
}

void StratifiedSample2D(float *samp, int nx, int ny, RNG &rng, bool jitter) {
    float dx = 1.f / nx, dy = 1.f / ny;
    for (int y = 0; y < ny; ++y)
        for (int x = 0; x < nx; ++x) {
            float jx = jitter ? rng.RandomFloat() : 0.5f;
            float jy = jitter ? rng.RandomFloat() : 0.5f;
            *samp++ = min((x + jx) * dx, OneMinusEpsilon);
            *samp++ = min((y + jy) * dy, OneMinusEpsilon);
        }
}

void LatinHypercube(float *samples, uint32_t nSamples, uint32_t nDim, RNG &rng) {
    // Generate LHS samples along diagonal
    float delta = 1.0f / nSamples;
    for (uint32_t i = 0; i < nSamples; ++i)
        for (uint32_t j = 0; j < nDim; ++j)
            samples[nDim * i + j] = min((i + (rng.RandomFloat())) * delta,
                                        OneMinusEpsilon);

    // Permute LHS samples in each dimension
    for (uint32_t i = 0; i < nDim; ++i) {
        for (uint32_t j = 0; j < nSamples; ++j) {
            uint32_t other = j + (rng.RandomUInt() % (nSamples - j));
            swap(samples[nDim * j + i], samples[nDim * other + i]);
        }
    }
}
