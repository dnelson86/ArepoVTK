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
		void CalculateMidpoints();
		void LimitCellDensities();
		
		// methods
		void DumpMesh();
		BBox WorldBound() const { return extent; }
    BBox VolumeBound() const { return extent; }

		// raster return
		bool TetraEdges(const int i, vector<Line> *edges);
    bool VoronoiEdges(const int i_face, vector<Line> *edges);
		
		// world geometry
    bool IntersectP(const Ray &r, float *t0, float *t1) const {
        return extent.IntersectP(r, t0, t1);
    }
		
		// mesh traversal
		void LocateEntryCell(const Ray &ray);
		void LocateEntryCellBrute(const Ray &ray);
		void VerifyPointInCell(int dp, Point &pos);
		
		int FindNearestGasParticle(Point &pt, double *mindist, int guess, int use_periodic);
		bool AdvanceRayOneCellNew(const Ray &ray, float *t0, float *t1, int previous_cell, 
															Spectrum &Lv, Spectrum &Tr);
		
		inline int getSphPID(int dp_id);
		
		// fluid data introspection
		inline double nnInterpScalar(int SphP_ID, int DP_ID, Vector &pt);
		inline void computeAuxVolumes(double *vol);
		float valMean(int valNum) { return valBounds[valNum*3+0]; }
		
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
    // rendering
		BBox extent;
    const TransferFunction *transferFunction;
		float viStepSize;
		
		// units, etc
		float valBounds[TF_NUM_VALS*3];     // min,max,mean for each non-derived quantity
		float unitConversions[TF_NUM_VALS]; // mult factor from code units to ArepoVTK "units"
		
		// mesh
    tessellation *T;
		tessellation AuxMesh;
		
		int *EdgeList;
		int *Nedges;
		int *NedgesOffset;
		
		// sunrise alternative connectivity construction
		
		vector<int>            primary_cells;   // maps dp_idx to the dp_idx of the primary cell
		vector<pair<int,int> > midpoint_idx;    // stores the starting midpoint and number of faces for the cell,
		                                       // indexed by dp_idx
		vector<Vector>         midpoints;       // stores the midpoints between the mesh points, indexed by midpoint_idx
		vector<int>            opposite_points; // stores the dp index of the cell opposite to the midpoint
};

#endif //ENABLE_AREPO
#endif //AREPO_RT_AREPO_H
