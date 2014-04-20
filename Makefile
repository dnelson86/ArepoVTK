#
# ArepoVTK
# Dylan Nelson
#

# user-configurable
# -----------------

EXECNAME = ArepoRT

#OPT += -DDEBUG
#OPT  += -DENABLE_OPENGL
#OPT  += -DENABLE_CUDA

# system
# -------

ARCH = $(shell uname)

OPTIMIZE = -Wall -g -m64 -O3 #-pg #enable profiler
INCL     = -I ./Arepo/

CC       = mpic++
LIBS     = -fopenmp -lm 

CFLAGS   = $(OPTIMIZE) -DH5_USE_16_API

ifeq ("$(SYSTYPE)","stampede_illustris")
CFLAGS += -I/home1/01637/mvogelsb/libs/gmp/include -I$(TACC_GSL_INC) -I$(TACC_HDF5_INC)
LIBS += -L/home1/01637/mvogelsb/libs/gmp/lib -L$(TACC_GSL_LIB) -L$(TACC_HDF5_LIB)
endif

ifeq ("$(SYSTYPE)","Gordon")
CFLAGS += -I/opt/gsl/gnu/include -I/opt/gnu/gmp/4.3.2/include -I/opt/hdf5/intel/mvapich2/ib/include
LIBS += -L/opt/gsl/gnu/lib -L/opt/gnu/gmp/4.3.2/lib -L/opt/hdf5/intel/mvapich2/ib/lib
endif

OBJS = ArepoRT.o camera.o fileio.o geometry.o integrator.o keyframe.o renderer.o sampler.o transfer.o transform.o util.o volume.o snapio.o
HEAD = ArepoRT.h camera.h fileio.h geometry.h integrator.h keyframe.h renderer.h sampler.h spectrum.h transfer.h transform.h util.h volume.h snapio.h
MISC_RM = frame.raw.txt frame.tga

# ENABLE_AREPO
OBJS += arepo.o arepoTree.o arepoInterp.o voronoi_3db.o
HEAD += arepo.h arepoTree.h
LIBS += -lgsl -lgslcblas -lgmp -lhdf5 -pthread -L ./Arepo/ -larepo #-lhwloc

OBJS := $(addprefix build/,$(OBJS))
INCL := $(addprefix src/,$(INCL))

$(EXECNAME): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(OPT) -o $(EXECNAME) $(LIBS)

clean:
	rm -f $(OBJS) $(EXECNAME) $(MISC_RM)

build/%.o: src/%.cpp
	$(CC) $(CFLAGS) $(OPT) -c $< -o $@
