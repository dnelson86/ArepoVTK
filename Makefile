#
# ArepoVTK
# Dylan Nelson
#

# user-configurable
# -----------------

EXECNAME = ArepoRT

#OPT += -DDEBUG          # enable verbose diagnostics and checks
#OPT  += -DENABLE_OPENGL # unused
#OPT  += -DENABLE_CUDA   # unused

# system
# -------

ARCH = $(shell uname)

OPTIMIZE = -Wall -g -m64 -O3 #-pg #enable profiler
INCL     = -I ./arepo/

CC       = mpicxx
LIBS     = -fopenmp -lm -L ./arepo/ 

CFLAGS   = $(OPTIMIZE) -DH5_USE_16_API

# module load gcc gsl fftw-serial hdf5-serial impi
CFLAGS += -I${GSL_HOME}/include -I${HDF5_HOME}/include -I./libpng/
LIBS += -L${GSL_HOME}/lib -L${HDF5_HOME}/lib -L./libpng/

OBJS = ArepoRT.o camera.o fileio.o fileio_img.o geometry.o integrator.o keyframe.o renderer.o sampler.o transfer.o transform.o util.o volume.o snapio.o
HEAD = ArepoRT.h camera.h fileio.h fileio_img.h geometry.h integrator.h keyframe.h renderer.h sampler.h spectrum.h transfer.h transform.h util.h volume.h snapio.h
MISC_RM = frame.raw.txt frame.tga

# ENABLE_AREPO
OBJS += arepo.o arepoTree.o arepoInterp.o voronoi_3db.o
HEAD += arepo.h arepoTree.h
LIBS += -lgsl -lgslcblas -lgmp -lhdf5 -pthread -larepo -lpng16 #-lhwloc

OBJS := $(addprefix build/,$(OBJS))
INCL := $(addprefix src/,$(INCL))

$(EXECNAME): libarepo.a $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(OPT) -o $(EXECNAME) $(LIBS)

libarepo.a:
	@cd arepo; make libarepo.a;

clean:
	@cd arepo; make clean;
	rm -f $(OBJS) $(EXECNAME) $(MISC_RM)

build/%.o: src/%.cpp
	$(CC) $(CFLAGS) $(OPT) -c $< -o $@
