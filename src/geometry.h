/*
 * geometry.h
 * dnelson
 */
 
#ifndef AREPO_RT_GEOMETRY_H
#define AREPO_RT_GEOMETRY_H

#include "voronoi_3db.h"

class Ray;

class Vector {
public:
    // construction
    Vector() {
			x = y = z = 0.0f;
		}
    Vector(double xx, double yy, double zz)
        : x(xx), y(yy), z(zz) {
    }
		Vector(float *xyz)
				: x(xyz[0]), y(xyz[1]), z(xyz[2]) {
		}
		Vector(double *xyz)
				: x(xyz[0]), y(xyz[1]), z(xyz[2]) {
		}
		explicit Vector(const Point &p);
		
		// operators
    Vector operator+(const Vector &v) const {
        return Vector(x + v.x, y + v.y, z + v.z);
    }
    Vector& operator+=(const Vector &v) {
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    Vector operator-(const Vector &v) const {
        return Vector(x - v.x, y - v.y, z - v.z);
    }
    Vector& operator-=(const Vector &v) {
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }
    Vector operator*(double f) const {
		  return Vector(f*x, f*y, f*z);
		}
    Vector &operator*=(double f) {
        x *= f; y *= f; z *= f;
        return *this;
    }
    Vector operator/(double f) const {
        double inv = 1.f / f;
        return Vector(x * inv, y * inv, z * inv);
    }
    Vector &operator/=(double f) {
        double inv = 1.f / f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    Vector operator-() const {
		  return Vector(-x, -y, -z);
		}
		
		void print(const string &preamble) const { cout << preamble << " x= " << setw(5) << x
		     << " y= " << setw(5) << y << " z= " << setw(5) << z << endl;
		}
		
    double operator[](int i) const { return (&x)[i]; }
    double &operator[](int i)      { return (&x)[i]; }
    inline double LengthSquared()   const { return x*x + y*y + z*z; }
    inline double Length()          const { return sqrt(LengthSquared()); }
		
    inline double PeriodicLengthSquared()   const {
			double xtmp, ytmp, ztmp;
			double dx = NGB_PERIODIC_LONG_X(x);
			double dy = NGB_PERIODIC_LONG_Y(y);
			double dz = NGB_PERIODIC_LONG_Z(z);
			return dx*dx + dy*dy + dz*dz;
		}
    inline double PeriodicLength()          const { return sqrt(PeriodicLengthSquared()); }

    bool operator==(const Vector &v) const { return x == v.x && y == v.y && z == v.z; }
    bool operator!=(const Vector &v) const { return x != v.x || y != v.y || z != v.z; }		
		
		// data
		double x, y, z;
};

class Point {
public:
    // construction
    Point() {
			x = y = z = 0.0;
		}
    Point(double xx, double yy, double zz)
        : x(xx), y(yy), z(zz) {
    }
		Point(float v[3])
				: x(v[0]), y(v[1]), z(v[2]) {
		}
		Point(double v[3])
				: x(v[0]), y(v[1]), z(v[2]) {
		}
		
		// operators
    Point operator+(const Vector &v) const {
        return Point(x + v.x, y + v.y, z + v.z);
    }
    Point &operator+=(const Vector &v) {
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    Vector operator-(const Point &p) const {
        return Vector(x - p.x, y - p.y, z - p.z);
    }
    Point operator-(const Vector &v) const {
        return Point(x - v.x, y - v.y, z - v.z);
    }
    Point &operator-=(const Vector &v) {
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }
    Point &operator+=(const Point &p) {
				x += p.x; y += p.y; z += p.z;
        return *this;
    }
    Point operator+(const Point &p) const {
        return Point(x + p.x, y + p.y, z + p.z);
    }
    Point operator* (double f) const {
        return Point(f*x, f*y, f*z);
    }
    Point &operator*=(double f) {
        x *= f; y *= f; z *= f;
        return *this;
    }
    Point operator/ (double f) const {
        double inv = 1.f/f;
        return Point(inv*x, inv*y, inv*z);
    }
    Point &operator/=(double f) {
        double inv = 1.f/f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
		
    double operator[](int i) const { return (&x)[i]; }
    double &operator[](int i) { return (&x)[i]; }
		
		void print(const string &preamble) const { cout << preamble << " x= " << setw(5) << x
		     << " y= " << setw(5) << y << " z= " << setw(5) << z << endl;
		}
		
    bool operator==(const Point &p) const { return x == p.x && y == p.y && z == p.z; }
    bool operator!=(const Point &p) const { return x != p.x || y != p.y || z != p.z; }
		
		// data
		double x, y, z;
};

class Ray {
public:
    Ray() : min_t(0.0f), max_t(INFINITY), index(-1), task(-1), prev_index(-1), depth(0), time(0.0f) { }
		
    Ray(const Point &origin, const Vector &direction, double start, double end = INFINITY, 
		    double t = 0.0f, int d = 0, int ind = -6, int tsk = -1, int pind = -1)
        : o(origin), d(direction), min_t(start), max_t(end), index(ind), 
				  task(tsk), prev_index(pind), depth(d), time(t) {
		}
		
    Point operator()(double t) const { return o + d * t; }
		void printRay(const string &preamble) { cout << preamble << " o.x= " << setw(5) << o.x
		     << " o.y= " << setw(5) << o.y << " o.z= " << setw(5) << o.z << " d.x= " << setw(5) << d.x 
				 << " d.y= " << setw(5) << d.y << " d.z= " << setw(5) << d.z << endl;
		}

    // data
    Point o;  // origin
    Vector d; // direction
    mutable double min_t, max_t; // parametric distance along ray
		mutable int index, task;    // index and task of the primary voronoi cell in which this ray is located
		mutable int prev_index;     // index of the primary voronoi cell where the ray was last
    mutable int depth; // counter of ViStepSize steps
		mutable int tetra; // DT index of the Delaunay tetra this ray is in (or was in last)
    double time;
};

class Line {
public:
    // construction
    Line(const Point &start, const Point&end)
        : p1(start), p2(end) { }

    // methods
    void Draw();

    // data
    Point p1;
    Point p2;

};

/* BBox */
class BBox {
public:
		// construction
    BBox() {
        pMin = Point( INFINITY,  INFINITY,  INFINITY);
        pMax = Point(-INFINITY, -INFINITY, -INFINITY);
    }
    BBox(const Point &p) : pMin(p), pMax(p) { }
    BBox(const Point &p1, const Point &p2) {
        pMin = Point(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z));
        pMax = Point(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z));
    }
		
		// inlines
    friend BBox Union(const BBox &b, const Point &p);
    friend BBox Union(const BBox &b, const BBox &b2);
    bool Overlaps(const BBox &b) const {
        bool x = (pMax.x >= b.pMin.x) && (pMin.x <= b.pMax.x);
        bool y = (pMax.y >= b.pMin.y) && (pMin.y <= b.pMax.y);
        bool z = (pMax.z >= b.pMin.z) && (pMin.z <= b.pMax.z);
        return (x && y && z);
    }
    bool Inside(const Point &pt) const {
        return (pt.x >= pMin.x && pt.x <= pMax.x &&
                pt.y >= pMin.y && pt.y <= pMax.y &&
                pt.z >= pMin.z && pt.z <= pMax.z);
    }
    void Expand(float delta) {
        pMin -= Vector(delta, delta, delta);
        pMax += Vector(delta, delta, delta);
    }
    float SurfaceArea() const {
        Vector d = pMax - pMin;
        return 2.0f * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    float Volume() const {
        Vector d = pMax - pMin;
        return d.x * d.y * d.z;
    }
    int MaximumExtent() const {
        Vector diag = pMax - pMin;
        if (diag.x > diag.y && diag.x > diag.z)
            return 0;
        else if (diag.y > diag.z)
            return 1;
        else
            return 2;
    }
    const Point &operator[](int i) const;
    Point &operator[](int i);
    Point Lerp(float tx, float ty, float tz) const {
        return Point(::Lerp(tx, pMin.x, pMax.x), ::Lerp(ty, pMin.y, pMax.y),
                     ::Lerp(tz, pMin.z, pMax.z));
    }
    Vector Offset(const Point &p) const {
        return Vector((p.x - pMin.x) / (pMax.x - pMin.x),
                      (p.y - pMin.y) / (pMax.y - pMin.y),
                      (p.z - pMin.z) / (pMax.z - pMin.z));
    }
		bool Edges(vector<Line> *edges);
    //void BoundingSphere(Point *c, float *rad) const;
    bool IntersectP(const Ray &ray, double *hitt0 = NULL, double *hitt1 = NULL) const;

		void print(const string &preamble) const {
				cout << preamble << " min.x= " << setw(5) << pMin.x << " min.y= " << setw(5) << pMin.y
				                 << " min.z= " << setw(5) << pMin.z << " max.x= " << setw(5) << pMax.x
												 << " max.y= " << setw(5) << pMax.y << " max.z= " << setw(5) << pMax.z << endl;
		}		
		
    bool operator==(const BBox &b) const { return b.pMin == pMin && b.pMax == pMax; }
    bool operator!=(const BBox &b) const { return b.pMin != pMin || b.pMax != pMax; }

    // data
    Point pMin, pMax;
};

// geometry inline functions
inline Vector::Vector(const Point &p)
{
		x = p.x;
		y = p.y;
		z = p.z;
}
inline Vector operator*(double f, const Vector &v) {
		return v*f;
}

// vector vector
inline double Dot(const Vector &v1, const Vector &v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
inline double AbsDot(const Vector &v1, const Vector &v2) {
    return fabsf(Dot(v1, v2));
}
inline Vector Cross(const Vector &v1, const Vector &v2) {
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
}

// vector
inline Vector Normalize(const Vector &v) {
		return v / v.Length();
}
inline void CoordinateSystem(const Vector &v1, Vector *v2, Vector *v3) {
    if (fabsf(v1.x) > fabsf(v1.y)) {
        double invLen = 1.0 / sqrtf(v1.x*v1.x + v1.z*v1.z);
        *v2 = Vector(-v1.z * invLen, 0.0f, v1.x * invLen);
    }
    else {
        double invLen = 1.0 / sqrtf(v1.y*v1.y + v1.z*v1.z);
        *v2 = Vector(0.0f, v1.z * invLen, -v1.y * invLen);
    }
    *v3 = Cross(v1, *v2);
}

// point
inline double Distance(const Point &p1, const Point &p2) {
    return (p1 - p2).Length();
}
inline double DistanceSquared(const Point &p1, const Point &p2) {
    return (p1 - p2).LengthSquared();
}
inline Point operator*(double f, const Point &p) {
    return p*f;
}

// BBox
inline const Point &BBox::operator[](int i) const {
    return (&pMin)[i];
}
inline Point &BBox::operator[](int i) {
    return (&pMin)[i];
}


#endif
