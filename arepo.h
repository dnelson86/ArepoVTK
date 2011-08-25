/*
 * arepo.h
 * dnelson
 */
 
#ifndef AREPO_RT_AREPO_H
#define AREPO_RT_AREPO_H
#ifdef ENABLE_AREPO

#include "geometry.h"
#include "transform.h"
#include "spectrum.h"
#include "volume.h"
#include "transfer.h"

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
class ArepoMesh {
public:
		// construction
		ArepoMesh(const TransferFunction *tf);
    ~ArepoMesh() {
        if (EdgeList)     delete EdgeList;
        if (Nedges)       delete Nedges;
        if (NedgesOffset) delete NedgesOffset;
    }

		// preprocessing
		void ComputeVoronoiEdges();
		void ComputeQuantityBounds();
		
		// methods
		void DumpMesh();
		BBox WorldBound() const { return extent; }
    BBox VolumeBound() const { return extent; }

		// raster return
		bool TetraEdges(const int i, vector<Line> *edges);
    bool VoronoiEdges(const int i_face, vector<Line> *edges);
		
    bool IntersectP(const Ray &r, float *t0, float *t1) const {
        return extent.IntersectP(r, t0, t1);
    }
		
		void LocateEntryCell(const Ray &ray, float *t0, float *t1);
		void LocateEntryCellBrute(const Ray &ray, float *t0, float *t1);
		int FindNearestGasParticle(Point &p, double *mindist);
		bool AdvanceRayOneCell(const Ray &ray, float *t0, float *t1, Spectrum &Lv, Spectrum &Tr);
		
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
		
private:
		BBox extent;
    const TransferFunction *transferFunction;
		float viStepSize;
		
		float densBounds[3]; // min,max,mean
		
    tessellation *T;
		
		int *EdgeList;
		int *Nedges;
		int *NedgesOffset;
};

#endif //ENABLE_AREPO
#endif //AREPO_RT_AREPO_H
