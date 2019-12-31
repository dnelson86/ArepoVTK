/*
 * arepoTree.h
 * dnelson
 */
 
#ifndef AREPO_RT_AREPOTREE_H
#define AREPO_RT_AREPOTREE_H

#include "transfer.h"

#if (NUM_THREADS > 1)
#include <omp.h>
#endif

// ArepoTree: expose the neighbor tree data structures and encapsulate tree related functions
class ArepoTree {
public:
  // construction
  ArepoTree(const TransferFunction *tf);
  ~ArepoTree();
  
  // methods
  BBox WorldBound() const { return extent; }
  BBox VolumeBound() const { return extent; }

  // world geometry
  bool IntersectP(const Ray &r, double *t0, double *t1) const {
    return extent.IntersectP(r, t0, t1);
  }
  
  // tree search traversal
  bool FindNeighborList(Point &pt, float hsml, int *numngb_int, vector<float> &vals);
  
  // sampling / interpolation
  bool AdvanceRayOneStep(const Ray &ray, double *t0, double *t1, Spectrum &Lv, Spectrum &Tr, int threadNum);
  
  // data
  
private:
  // rendering
  BBox extent;
  const TransferFunction *transferFunction;
  
  //vector<int *> varNGBLists;
  //int *varNGBList;
};

#endif //AREPO_RT_AREPOTREE_H
