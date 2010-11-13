
OPTFLAG = -O3
CXX = mpic++
FC = mpif77
CC = mpic++
FFLAGS += $(OPTFLAG)
CXXFLAGS += -Wall -DBLACS_ALL -Iinclude -Wno-write-strings $(OPTFLAG)
LDLIBS += -lscalapack-openmpi -lblacsCinit-openmpi -lblacs-openmpi
#  -llapack -lblas -latlas

PSTEGR_SRC := disnan.o   dlar1v.o  dlarrb2.o  dlarrc.o   dlarrd.o    \
              dlarre2.o  dlarrf.o  dlarrv2.o  dstegr2a.o  dstegr2.o    \
              dlar1va.o  dlarra.o  dlarrb.o   dlarrd2.o  dlarre2a.o    \
              dlarrf2.o  dlarrk.o  dlarrv.o   dstegr2b.o  
             
PSTEGR_DIR := pdsyevr/pstegr
PSTEGR := $(addprefix $(PSTEGR_DIR)/,$(PSTEGR_SRC))

all : lapeigs graph_eigs

.PHONY : clean test

lapeigs.o :  scalapack_symmetric_eigen.hpp scalapack_symmetric_eigen.cc

lapeigs :  lapeigs.o mpiutil.o  pdsyevr/pdsyevr.o pdsyevr/pdsyevr_tri.o $(PSTEGR)

graph_eigs.o : graph_eigs_opts.hpp triplet_graph.hpp scalapack_symmetric_eigen.cc

graph_eigs : graph_eigs.o mpiutil.o pdsyevr/pdsyevr.o pdsyevr/pdsyevr_tri.o $(PSTEGR)

clean:
	rm -rf graph_eigs graph_eigs.o lapeigs lapeigs.o mpiutil.o  pdsyevr/pdsyevr.o pdsyevr/pdsyevr_tri.o $(PSTEGR)

test: lapeigs
	./lapeigs test/tapir.smat test.evals N
