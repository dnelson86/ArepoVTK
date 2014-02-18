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
##CPPFLAGS = $(OPT) $(INCL)

OBJS = ArepoRT.o camera.o fileio.o geometry.o integrator.o keyframe.o renderer.o sampler.o transfer.o transform.o util.o volume.o snapio.o
HEAD = ArepoRT.h camera.h fileio.h geometry.h integrator.h keyframe.h renderer.h sampler.h spectrum.h transfer.h transform.h util.h volume.h snapio.h
MISC_RM = frame.raw.txt frame.tga

# ENABLE_AREPO
OBJS += arepo.o arepoTree.o arepoInterp.o voronoi_3db.o
HEAD += arepo.h arepoTree.h
LIBS += -lgsl -lgslcblas -lgmp -lhdf5 -lhwloc -pthread -L ./Arepo/ -larepo

OBJS := $(addprefix build/,$(OBJS))
INCL := $(addprefix src/,$(INCL))

$(EXECNAME): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(OPT) -o $(EXECNAME) $(LIBS)

clean:
	rm -f $(OBJS) $(EXECNAME) $(MISC_RM)

build/%.o: src/%.cpp
	$(CC) $(CFLAGS) $(OPT) -c $< -o $@
