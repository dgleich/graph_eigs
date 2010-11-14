
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

.PHONY : clean test clean_test

lapeigs.o :  scalapack_symmetric_eigen.hpp scalapack_symmetric_eigen.cc

lapeigs :  lapeigs.o mpiutil.o  pdsyevr/pdsyevr.o pdsyevr/pdsyevr_tri.o $(PSTEGR)

graph_eigs.o : graph_eigs_opts.hpp triplet_graph.hpp scalapack_symmetric_eigen.cc

graph_eigs : graph_eigs.o mpiutil.o pdsyevr/pdsyevr.o pdsyevr/pdsyevr_tri.o $(PSTEGR)

clean:
	rm -rf graph_eigs graph_eigs.o lapeigs lapeigs.o mpiutil.o  pdsyevr/pdsyevr.o pdsyevr/pdsyevr_tri.o $(PSTEGR)
	
clean_test:
	rm -rf test/tapir.smat.*
	rm -rf test/polblogs-sym-cc.smat.*
	
test_tapir:	
	./graph_eigs test/tapir.smat -a -r -p > /dev/null
	python graph_eigs.py test/tapir.smat -a --check-eigs test/tapir.smat.adjacency.eigs --check-ipar test/tapir.smat.adjacency.ipar
	./graph_eigs test/tapir.smat -l -r -p > /dev/null
	python graph_eigs.py test/tapir.smat -l --check-eigs test/tapir.smat.laplacian.eigs --check-ipar test/tapir.smat.laplacian.ipar
	./graph_eigs test/tapir.smat -n -r -p > /dev/null
	python graph_eigs.py test/tapir.smat -n --check-eigs test/tapir.smat.normalized.eigs --check-ipar test/tapir.smat.normalized.ipar
	./graph_eigs test/tapir.smat -m -r -p > /dev/null
	python graph_eigs.py test/tapir.smat -m --check-eigs test/tapir.smat.modularity.eigs --check-ipar test/tapir.smat.modularity.ipar
	
test_polblogs:
	./graph_eigs test/polblogs-sym-cc.smat -a -r -p > /dev/null
	python graph_eigs.py test/polblogs-sym-cc.smat -a --check-eigs test/polblogs-sym-cc.smat.adjacency.eigs --check-ipar test/polblogs-sym-cc.smat.adjacency.ipar
	./graph_eigs test/polblogs-sym-cc.smat -l -r -p > /dev/null
	python graph_eigs.py test/polblogs-sym-cc.smat -l --check-eigs test/polblogs-sym-cc.smat.laplacian.eigs --check-ipar test/polblogs-sym-cc.smat.laplacian.ipar
	./graph_eigs test/polblogs-sym-cc.smat -n -r -p > /dev/null
	python graph_eigs.py test/polblogs-sym-cc.smat -n --check-eigs test/polblogs-sym-cc.smat.normalized.eigs --check-ipar test/polblogs-sym-cc.smat.normalized.ipar
	./graph_eigs test/polblogs-sym-cc.smat -m -r -p > /dev/null
	python graph_eigs.py test/polblogs-sym-cc.smat -m --check-eigs test/polblogs-sym-cc.smat.modularity.eigs --check-ipar test/polblogs-sym-cc.smat.modularity.ipar

test: graph_eigs clean_test test_tapir test_polblogs

	
