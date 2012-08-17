/*
 * volume.h
 * dnelson
 */
 
#ifndef AREPO_RT_VOLUME_H
#define AREPO_RT_VOLUME_H

#include "ArepoRT.h"
#include "transform.h"

class VolumeRegion {
public:
    // construction
    virtual ~VolumeRegion();
		
		// pure virtual methods
    virtual BBox WorldBound() const = 0;
    virtual bool IntersectP(const Ray &ray, float *t0, float *t1) const = 0;
    virtual Spectrum sigma_a(const Point &, const Vector &, float time) const = 0;
    virtual Spectrum sigma_s(const Point &, const Vector &, float time) const = 0;
    virtual Spectrum Lve(const Point &, const Vector &, float time) const = 0;
    virtual Spectrum sigma_t(const Point &p, const Vector &wo, float time) const = 0;
    virtual Spectrum tau(const Ray &ray, float step = 1.f, float offset = 0.5) const = 0;
};

#include "arepo.h"

// Scene
class Scene {
public:
    // construction
    Scene(VolumeRegion *vr, ArepoMesh* am);
    ~Scene();
		
    //bool Intersect(const Ray &ray, Intersection *isect) const {
    //    bool hit = volumeRegion->Intersect(ray, isect);
    //    return hit;
    //}
    bool IntersectP(const Ray &ray, float *t0, float *t1) const {
        bool hit = volumeRegion->IntersectP(ray, t0, t1);
        return hit;
    }
    const BBox &WorldBound() const;

    // data
    VolumeRegion *volumeRegion;
		ArepoMesh *arepoMesh;
    BBox bound;
};

class DensityRegion : public VolumeRegion {
public:
		// construction
    DensityRegion(const Spectrum &sa, const Spectrum &ss, 
                  const Spectrum &emit, const Transform &VolumeToWorld)
        : sig_a(sa), sig_s(ss), le(emit), WorldToVolume(Inverse(VolumeToWorld)) { }
					
		// virtuals
    virtual float Density(const Point &Pobj) const = 0;
		
		// defined as simple scalings with density
    Spectrum sigma_a(const Point &p, const Vector &, float) const {
        return Density(WorldToVolume(p)) * sig_a;
    }
    Spectrum sigma_s(const Point &p, const Vector &, float) const {
        return Density(WorldToVolume(p)) * sig_s;
    }
    Spectrum sigma_t(const Point &p, const Vector &, float) const {
        return Density(WorldToVolume(p)) * (sig_a + sig_s);
    }
    Spectrum Lve(const Point &p, const Vector &, float) const {
        return Density(WorldToVolume(p)) * le;
    }
    Spectrum tau(const Ray &r, float stepSize, float offset) const;
		
protected:
    // data
    Spectrum sig_a, sig_s, le;
    Transform WorldToVolume;
};

class VolumeGridDensity : public DensityRegion {
public:
    // construction
    VolumeGridDensity(const Spectrum &sa, const Spectrum &ss, 
            const Spectrum &emit, const BBox &e, const Transform &v2w,
            int x, int y, int z, const float *d)
        : DensityRegion(sa, ss, emit, v2w), nx(x), ny(y), nz(z), extent(e)
    {
				IF_DEBUG(cout << "VolumeGridDensity(" << nx << ", " << ny << ", " << nz << ") constructor." << endl);

        density = new float[nx*ny*nz];
        memcpy(density, d, nx*ny*nz*sizeof(float));
    }
    ~VolumeGridDensity() { delete[] density; }
		
		// methods
    BBox WorldBound() const { return Inverse(WorldToVolume)(extent); }
    BBox VolumeBound() const { return extent; }
    bool IntersectP(const Ray &r, float *t0, float *t1) const {
				//IF_DEBUG(*r->printRay("DR IntersectP W "));
        Ray ray = WorldToVolume(r);
        return extent.IntersectP(ray, t0, t1);
    }
    float Density(const Point &Pobj) const;
    float D(int x, int y, int z) const {
        x = Clamp(x, 0, nx-1);
        y = Clamp(y, 0, ny-1);
        z = Clamp(z, 0, nz-1);
        return density[z*nx*ny + y*nx + x];
    }
private:
    // data
    float *density;
    const int nx, ny, nz;
    const BBox extent;
};

VolumeGridDensity *CreateGridVolumeRegion(const Transform &volume2world, const string &filename);

#endif
