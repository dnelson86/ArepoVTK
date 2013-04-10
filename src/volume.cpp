/*
 * volume.cpp
 * dnelson
 */

#include "volume.h"
#include "spectrum.h"

// Scene
Scene::Scene(VolumeRegion *vr, ArepoMesh *am)
{
		IF_DEBUG(cout << "Scene() constructor." << endl);
		
    volumeRegion = vr;
		arepoMesh    = am;
		
    if (volumeRegion) bound = volumeRegion->WorldBound();
    if (arepoMesh) bound = Union(bound, arepoMesh->WorldBound());
}

Scene::~Scene()
{
    delete volumeRegion;
		delete arepoMesh;
}

const BBox &Scene::WorldBound() const
{
    return bound;
}


VolumeRegion::~VolumeRegion()
{
}

Spectrum DensityRegion::tau(const Ray &r, float stepSize, float u) const
{
    double t0, t1;
    float length = r.d.Length();
    if (length == 0.f) return 0.f;
    Ray rn(r.o, r.d / length, r.min_t * length, r.max_t * length, r.time);
    if (!IntersectP(rn, &t0, &t1)) return 0.;
    Spectrum tau(0.0);
    t0 += u * stepSize;
    while (t0 < t1) {
        tau += sigma_t(rn(t0), -rn.d, r.time);
        t0 += stepSize;
    }
    return tau * stepSize;
}

float VolumeGridDensity::Density(const Point &Pobj) const
{
    if (!extent.Inside(Pobj)) return 0;
		
    // compute voxel coordinates and offsets for Pobj
    Vector vox = extent.Offset(Pobj);
    vox.x = vox.x * nx - 0.5f;
    vox.y = vox.y * ny - 0.5f;
    vox.z = vox.z * nz - 0.5f;
		
    int vx = (int)(vox.x);
		int vy = (int)(vox.y);
		int vz = (int)(vox.z);
		
    float dx = vox.x - vx;
		float dy = vox.y - vy;
		float dz = vox.z - vz;

		IF_DEBUG(cout << " Density: vox.x = " << vox.x << " vox.y = " << vox.y << " vox.z = " << vox.z 
									<< " vx = " << vx << " vy = " << vy << " vz = " << vz << endl);

    // Trilinearly interpolate density values to compute local density
    float d00 = Lerp(dx, D(vx, vy, vz),     D(vx+1, vy, vz));
    float d10 = Lerp(dx, D(vx, vy+1, vz),   D(vx+1, vy+1, vz));
    float d01 = Lerp(dx, D(vx, vy, vz+1),   D(vx+1, vy, vz+1));
    float d11 = Lerp(dx, D(vx, vy+1, vz+1), D(vx+1, vy+1, vz+1));
    float d0 = Lerp(dy, d00, d10);
    float d1 = Lerp(dy, d01, d11);

		IF_DEBUG(cout << " Density: d00 = " << d00 << " d10 = " << d10 << " d01 = " << d01 << " d11 = "
									<< d11 << " d0 = " << d0 << " d1 = " << d1 << endl);
    return Lerp(dz, d0, d1);
}

VolumeGridDensity *CreateGridVolumeRegion(const Transform &volume2world, const string &filename)
{
    // Initialize common volume region parameters
    Spectrum sigma_a = 0.0;
    Spectrum sigma_s = 0.0;
    Spectrum Le = Spectrum::FromNamed("red");
    Point p0 = Point(0,0,0);
    Point p1 = Point(1,1,1);
		
    int ni,nx,ny,nz;
    vector<float> data;

	  if (!(ni = parseSceneFile(filename, nx, ny, nz, &data))) {
	    cout << "ERROR: Couldn't open scene file: " << filename << endl;
		  return 0;
	  }
		
		cout << "Read [" << ni << "] points from scene file: " << filename << endl;

    if (data.empty()) {
        cout << "No \"density\" values provided for volume grid?" << endl;
        return NULL;
    }
    if (ni != nx*ny*nz) {
        cout << "VolumeGridDensity has " << ni << " density values but nx*ny*nz = "
				     << nx*ny*nz << endl;
        return NULL;
    }

#if defined(DEBUG) && 0
    for (int i=0; i < data.size(); i++)
				cout << " CGVR data[" << i << "] = " << data[i] << endl;
#endif
		
    return new VolumeGridDensity(sigma_a, sigma_s, Le, BBox(p0, p1),
        volume2world, nx, ny, nz, (float*)&data[0]);
}

