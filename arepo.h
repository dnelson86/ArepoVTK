/*
 * arepo.h
 * dnelson
 */
 
#ifndef AREPO_RT_AREPO_H
#define AREPO_RT_AREPO_H
#ifdef ENABLE_AREPO

// Arepo overrides
// #define mpi_printf(expr) ((void)0)

#include "geometry.h"
#include "transform.h"
#include "spectrum.h"
#include "volume.h"

#include "mpi.h" //as in Sunrise, need to include outside the C-linkage block
#include "gmp.h"

extern "C" {
#include "arepoconfig.h"
#include "proto.h"   //allvars.h (mpi,gsl)
#include "voronoi.h" //gmp
}

// Arepo: main interface with Arepo to load a snapshot, create data structures, and return
class Arepo
{
public:
    // construction
    Arepo(const string &sfn, const string &pfn) : snapFilename(sfn), paramFilename(pfn) { }

    // methods
    void Init(int*, char***);
		void Cleanup();
    bool LoadSnapshot();

private:
		// data
    string snapFilename;
		string paramFilename;
};

// ArepoMesh: expose the Voronoi data structures and encapsulate mesh related functions
class ArepoMesh : public VolumeRegion {
public:
		// construction
		ArepoMesh(const Spectrum &sa, const Spectrum &ss, const Spectrum &emit, const Transform &VolumeToWorld);
    ~ArepoMesh() {
        if (EdgeList)     delete EdgeList;
        if (Nedges)       delete Nedges;
        if (NedgesOffset) delete NedgesOffset;
    }

		// preprocessing
		void ComputeVoronoiEdges();
		
		// methods
		void DumpMesh();
    BBox WorldBound() const { return Inverse(WorldToVolume)(extent); }
    BBox VolumeBound() const { return extent; }

		// raster return
		bool TetraEdges(const int i, vector<Line> *edges);
    bool VoronoiEdges(const int i_face, vector<Line> *edges);
		
    bool IntersectP(const Ray &r, float *t0, float *t1) const {
				//IF_DEBUG(*r->printRay("AM IntersectP W "));
        Ray ray = WorldToVolume(r);
        return extent.IntersectP(ray, t0, t1);
    }
		
		void LocateEntryCell(const Ray &ray, float *t0, float *t1);
		bool AdvanceRayOneCell(const Ray &ray, float *t0, float *t1, Spectrum &Lv, Spectrum &Tr);
		
		//TODO: change these from taking inputs=Point p, which would require we go into the mesh and calculate the
		//      hydro quantites at that point (extremely expensive). Rather, this is already accessible inside the
		//      integration stage, so pass these quantities in and have e.g. Lve encapsulate the transfer functions.
    Spectrum sigma_a(const Point &p, const Vector &, float) const {    }
    Spectrum sigma_s(const Point &p, const Vector &, float) const {    }
    Spectrum sigma_t(const Point &p, const Vector &, float) const {    }
    Spectrum Lve(const Point &p, const Vector &, float) const {    }
    Spectrum tau(const Ray &r, float stepSize, float offset) const {   }
		
		// data
		int Ndp;           // number of delaunay points
		int Ndt;           // number of delaunay tetrahedra
		int Nvf;           // number of voronoi faces		
		
		point *DP;         // delaunay points
		tetra *DT;         // delaunay tetrahedra
		tetra_center *DTC; // circumcenters of delaunay tetrahedra
		char *DTF;         // tetra faces
		face *VF;          // voronoi faces
    //connection *DC;    // voronoi connections		
		
		Transform WorldToVolume;
		
private:
		BBox extent;
    Spectrum sig_a, sig_s, le;
		
    tessellation *T;
		
		int *EdgeList;
		int *Nedges;
		int *NedgesOffset;
};

ArepoMesh *CreateArepoMesh(const Transform &volume2world);

#endif //ENABLE_AREPO
#endif //AREPO_RT_AREPO_H
