#!/bin/bash            # this line only there to enable syntax highlighting in this file

##################################################
#  Enable/Disable compile-time options as needed #
##################################################

#DISABLE_MEMORY_MANAGER
SUNRISE # disables Inf terminate in get_tetra
PERIODIC
VORONOI
CHUNKING
VORONOI_DYNAMIC_UPDATE # alters SphP
NO_FANCY_MPI_CONSTRUCTS
NO_MPI_IN_PLACE
NO_ISEND_IRECV_IN_DOMAIN
FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
HAVE_HDF5
VORONOI_MESHOUTPUT # for write_voronoi_mesh
#MESHRELAX_DENSITY_IN_INPUT # ICs (not for normal snapshots)
#INPUT_IN_DOUBLEPRECISION # DP ICs or DP snaps

# for ArepoVTK only (not Arepo projection)
#NUM_THREADS=4
#SPECIAL_BOUNDARY # for test.spoon, doesn't hurt otherwise (allows negative IDs)
#AREPOVTK_MINMEM_TREEONLY # NOT IMPLEMENTED FULLY # P/SphP only contain: Pos,Utherm,Density (all that is loaded)

# illustris.fof0 (for Arepo projection, not ArepoVTK)
VORONOI_NEW_IMAGE
COOLING
USE_SFR
METALS
