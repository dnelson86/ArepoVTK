#
# ArepoVTK
# Dylan Nelson
#

# user-configurable
# -----------------

EXECNAME = ArepoRT

#OPT += -DDEBUG
OPT  += -DENABLE_AREPO
#OPT  += -DENABLE_OPENGL
#OPT  += -DENABLE_CUDA

# system
# -------

ARCH = $(shell uname)

OPTIMIZE = -Wall -g -m64 -O3
INCL     = -I ./Arepo/

CC       = mpic++
LIBS     = -lm 

CFLAGS   = $(OPTIMIZE)
CPPFLAGS = $(OPT) $(INCL)

OBJS = ArepoRT.o camera.o fileio.o geometry.o integrator.o renderer.o sampler.o transfer.o transform.o util.o volume.o
HEAD = ArepoRT.h camera.h fileio.h geometry.h integrator.h renderer.h sampler.h spectrum.h transfer.h transform.h util.h volume.h
MISC_RM = frame.raw.txt frame.tga

ifeq (ENABLE_AREPO,$(findstring ENABLE_AREPO,$(OPT)))
  OBJS += arepo.o voronoi_3db.o
  HEAD += arepo.h
  LIBS += -lgsl -lgslcblas -lgmp -lhdf5 -pthread -L ./Arepo/ -larepo
endif

OBJS := $(addprefix build/,$(OBJS))
INCL := $(addprefix src/,$(INCL))

$(EXECNAME): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(OPT) -o $(EXECNAME) $(LIBS)

clean:
	rm -f $(OBJS) $(EXECNAME) $(MISC_RM)

build/%.o: src/%.cpp
	$(CC) $(CFLAGS) $(OPT) -c $< -o $@
