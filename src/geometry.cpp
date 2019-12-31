/*
 * geometry.cpp
 * dnelson
 */

#include "ArepoRT.h" 
#include "geometry.h"

//BBox

bool BBox::Edges(vector<Line> *edges)
{
  //edges->push_back(Line(Point(pMin.x,pMin.y,pMin.z),Point(pMax.x,pMax.y,pMax.z))); //LL->UR diag debug

  edges->push_back(Line(Point(pMin.x,pMin.y,pMin.z),Point(pMin.x,pMax.y,pMin.z))); //1 - FF Left
  edges->push_back(Line(Point(pMin.x,pMax.y,pMin.z),Point(pMax.x,pMax.y,pMin.z))); //2 - FF Top
  edges->push_back(Line(Point(pMax.x,pMax.y,pMin.z),Point(pMax.x,pMin.y,pMin.z))); //3 - FF Right
  edges->push_back(Line(Point(pMax.x,pMin.y,pMin.z),Point(pMin.x,pMin.y,pMin.z))); //4 - FF Bottom
  edges->push_back(Line(Point(pMin.x,pMax.y,pMin.z),Point(pMin.x,pMax.y,pMax.z))); //5 - connect 1/2+12/9
  edges->push_back(Line(Point(pMax.x,pMax.y,pMin.z),Point(pMax.x,pMax.y,pMax.z))); //6 - connect 2/3+9/10
  edges->push_back(Line(Point(pMax.x,pMin.y,pMin.z),Point(pMax.x,pMin.y,pMax.z))); //7 - connect 3/4+10/11
  edges->push_back(Line(Point(pMin.x,pMin.y,pMin.z),Point(pMin.x,pMin.y,pMax.z))); //8 - connect 4/1+11/12
  edges->push_back(Line(Point(pMin.x,pMin.y,pMax.z),Point(pMin.x,pMax.y,pMax.z))); //9 - RF Left
  edges->push_back(Line(Point(pMin.x,pMax.y,pMax.z),Point(pMax.x,pMax.y,pMax.z))); //10 - RF Top
  edges->push_back(Line(Point(pMax.x,pMax.y,pMax.z),Point(pMax.x,pMin.y,pMax.z))); //11 - RF Right
  edges->push_back(Line(Point(pMax.x,pMin.y,pMax.z),Point(pMin.x,pMin.y,pMax.z))); //12 - RF Bottom
  
  return true;
}

bool BBox::IntersectP(const Ray &ray, double *hitt0,
                      double *hitt1) const
{
  double t0 = ray.min_t, t1 = ray.max_t;
  
  IF_DEBUG(cout << "BBox:IntersectP(t0 = " << t0 << " t1 = " << t1 << ") ray o.x=" << ray.o.x
                << " o.y=" << ray.o.y << " o.z=" << ray.o.z << " d.x=" << ray.d.x 
                << " d.y=" << ray.d.y << " d.z=" << ray.d.z << endl);
  
  for (int i = 0; i < 3; ++i) {
    // Update interval for _i_th bounding box slab
    double invRayDir = 1.0 / ray.d[i];
    double tNear = (pMin[i] - ray.o[i]) * invRayDir;
    double tFar  = (pMax[i] - ray.o[i]) * invRayDir;

    // Update parametric interval from slab intersection $t$s
    if (tNear > tFar) swap(tNear, tFar);
    t0 = tNear > t0 ? tNear : t0;
    t1 = tFar  < t1 ? tFar  : t1;
    IF_DEBUG(cout << " i[" << i << "] t0 = " << t0 << " t1 = " << t1 << endl);
    if (t0 > t1) return false;
  }
  if (hitt0) *hitt0 = t0;
  if (hitt1) *hitt1 = t1;
  return true;
}

BBox Union(const BBox &b, const Point &p) {
  BBox ret = b;
  ret.pMin.x = min(b.pMin.x, p.x);
  ret.pMin.y = min(b.pMin.y, p.y);
  ret.pMin.z = min(b.pMin.z, p.z);
  ret.pMax.x = max(b.pMax.x, p.x);
  ret.pMax.y = max(b.pMax.y, p.y);
  ret.pMax.z = max(b.pMax.z, p.z);
  return ret;
}


BBox Union(const BBox &b, const BBox &b2) {
  BBox ret;
  ret.pMin.x = min(b.pMin.x, b2.pMin.x);
  ret.pMin.y = min(b.pMin.y, b2.pMin.y);
  ret.pMin.z = min(b.pMin.z, b2.pMin.z);
  ret.pMax.x = max(b.pMax.x, b2.pMax.x);
  ret.pMax.y = max(b.pMax.y, b2.pMax.y);
  ret.pMax.z = max(b.pMax.z, b2.pMax.z);
  return ret;
}
