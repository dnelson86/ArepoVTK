/*
 * voronoi_3db.h
 * dnelson
 */
 
#ifndef VORONOI_3DB_H
#define VORONOI_3DB_H
 
#include "ArepoRT.h"
// otherwise this file doesn't seem to know about the preprocessor definitions from Config.sh
#include "../Arepo/build/arepoconfig.h"

#ifndef SPECIAL_BOUNDARY
#ifndef LONGIDS
typedef unsigned int MyIDType;
#define MPI_MYIDTYPE MPI_UNSIGNED
#else
typedef unsigned long long MyIDType;
#define MPI_MYIDTYPE MPI_UNSIGNED_LONG_LONG
#endif
#else
#ifndef LONGIDS
typedef int MyIDType;
#define MPI_MYIDTYPE MPI_INT
#else
typedef long long MyIDType;
#define MPI_MYIDTYPE MPI_LONG_LONG_INT
#endif
#endif

#ifndef DOUBLEPRECISION         // default is single-precision
typedef float MyFloat;
typedef float MyDouble;
#define MPI_MYFLOAT MPI_FLOAT
#define MPI_MYDOUBLE MPI_FLOAT
#else
#if (DOUBLEPRECISION == 2)      // mixed precision
typedef float MyFloat;
typedef double MyDouble;
#define MPI_MYFLOAT MPI_FLOAT
#define MPI_MYDOUBLE MPI_DOUBLE
#else // everything double-precision
typedef double MyFloat;
typedef double MyDouble;
#define MPI_MYFLOAT MPI_DOUBLE
#define MPI_MYDOUBLE MPI_DOUBLE
#endif
#endif

#define MAXGRADIENTS 5

#include "mpi.h" //as in Sunrise, need to include outside the C-linkage block
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <assert.h>
#include "gmp.h"

extern "C" {
#include "../Arepo/build/arepoconfig.h"
#include "../Arepo/src/mesh.h"
#include "../Arepo/src/voronoi.h" //gmp
}

#include "allvars.h"

void init_clear_auxmesh(tessellation * T);
int insert_point_new(tessellation * T, int pp, int ttstart);
void make_an_edge_split_new(tessellation * T, int tt0, int edge_nr, int count, int pp, int *ttlist);
void compute_auxmesh_volumes(tessellation *T, double *vol);
void derefine_refine_process_edge_new(tessellation * T, double *vol, int tt, int nr, unsigned char *visited_edges);

// for DC connectivity
int find_next_cell_DC(tessellation *T, int cell, double p0[3], double dir[3], int previous, double *length);

// for delaunay based NNI
bool calc_circumcenter(tessellation *T, point *p0, int dp1, int dp2, int dp3, double *cp);

#endif // VORONOI_3DB_H
