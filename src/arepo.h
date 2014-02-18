/*
 * arepo.h
 * dnelson
 */
 
#ifndef AREPO_RT_AREPO_H
#define AREPO_RT_AREPO_H

#ifdef DEUBG
#define VERBOSE
#endif

#include "transfer.h"
#include "voronoi_3db.h"

#if (NUM_THREADS > 1)
#include <omp.h>
#endif

void addValsContribution( vector<float> &vals, int SphP_ind, double weight );

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
    ~ArepoMesh();
		
		// init for particular interp methods
		void setupAuxMeshes();
		void precomputeTetraGrads();

		// preprocessing
		int ComputeVoronoiEdges();
		void ComputeQuantityBounds();
		void LimitCellDensities();
		
		// methods
		void DumpMesh();
		void OutputMesh();
		BBox WorldBound() const { return extent; }
    BBox VolumeBound() const { return extent; }

		// raster return
		bool TetraEdges(const int i, vector<Line> *edges);
    bool VoronoiEdges(const int i_face, vector<Line> *edges);
		
		// world geometry
    bool IntersectP(const Ray &r, double *t0, double *t1) const {
        return extent.IntersectP(r, t0, t1);
    }
		
		// mesh traversal
		void LocateEntryCell(const Ray &ray, int *prevEntryCell);
		void LocateEntryCellBrute(const Ray &ray);
		void VerifyPointInCell(int sphInd, Point &pos);
		
		void LocateEntryTetra(const Ray &ray, int *prevEntryTetra);
		
		int FindNearestGasParticle(Point &pt, int guess, double *mindist);
		bool AdvanceRayOneCellNew(const Ray &ray, double *t0, double *t1, 
															Spectrum &Lv, Spectrum &Tr, int threadNum);
		
		inline int getSphPID(int dpInd);
		void locateCurrentTetra(const Ray& ray, Vector &pt);
		void checkCurCellTF(bool *addFlag, int sphInd, vector<float> &vals);
		
		// fluid data introspection
		float calcNeighborHSML(int sphInd, Point &pt);
		int subSampleCell(const Ray &ray, Point &pt, vector<float> &vals, int threadNum);
		float valMean(int valNum) { return valBounds[valNum*3+0]; }
		
		// NNI_WATSON_SAMBRIDGE
		inline bool needTet(int tt, point *pp, int *node_inds, int *nTet);
		void addTet(int tt, point *pp, int *node_inds, int *tet_inds, int *nNode, int *nTet);
		double ccVolume(double *ci, double *cj, double *ck, double *ct);
		
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
		
		// units, etc
		float valBounds[TF_NUM_VALS*3];     // min,max,mean for each non-derived quantity
		float unitConversions[TF_NUM_VALS]; // mult factor from code units to ArepoVTK "units"
		
		// mesh
		tessellation *T;
		
		// for particular interpolation methods
		tessellation *AuxMeshes;
		float *DT_grad;
		float *DP_vols;
		
		// for drawing voronoi faces and edges
		vector<int> vertexList;
		vector<int> numVertices;
		vector<int> vertexOffset;
};

#endif //AREPO_RT_AREPO_H
