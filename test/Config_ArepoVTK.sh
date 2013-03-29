#!/bin/bash            # this line only there to enable syntax highlighting in this file

##################################################
#  Enable/Disable compile-time options as needed #
##################################################

#DISABLE_MEMORY_MANAGER
SUNRISE # disables Inf terminate in get_tetra
PERIODIC
VORONOI
#MESHRELAX_DENSITY_IN_INPUT (not for normal snapshots)
CHUNKING
VORONOI_DYNAMIC_UPDATE # alters SphP
NO_FANCY_MPI_CONSTRUCTS
NO_MPI_IN_PLACE
NO_ISEND_IRECV_IN_DOMAIN
FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
HAVE_HDF5
#VORONOI_MESHOUTPUT # for write_voronoi_mesh
#DOUBLEPRECISION=1 # only for test.diego right now
#INPUT_IN_DOUBLEPRECISION # only for test.diego right now
