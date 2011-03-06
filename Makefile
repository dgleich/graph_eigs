
OPTFLAG = -O3
CXX = mpic++
FC = mpif77
CC = mpic++
FFLAGS += $(OPTFLAG) -funroll-loops
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

graph_eigs : graph_eigs.o mpiutil.o pdsyevr/pdsyevr.o pdsyevr/pdsyevr_tri.o pdsyev_tri.o $(PSTEGR) 

clean:
	rm -rf graph_eigs graph_eigs.o lapeigs lapeigs.o mpiutil.o  pdsyevr/pdsyevr.o pdsyevr/pdsyevr_tri.o $(PSTEGR)
	
clean_test:
	rm -rf $(all_tests_clean)
define TEST_template 
test_$(1) : graph_eigs
	rm -rf test/$(1).smat.*
	$(2) ./graph_eigs test/$(1).smat -a -r -p > /dev/null
	python graph_eigs.py test/$(1).smat -a --check-eigs test/$(1).smat.adjacency.eigs --check-ipar test/$(1).smat.adjacency.ipar --check-resids test/$(1).smat.adjacency.resids
	$(2) ./graph_eigs test/$(1).smat -l -r -p --commute-all > /dev/null
	python graph_eigs.py test/$(1).smat -l \
	    --check-eigs test/$(1).smat.laplacian.eigs \
	    --check-ipar test/$(1).smat.laplacian.ipar \
	    --check-resids test/$(1).smat.laplacian.resids \
	    --check-commute-all test/$(1).smat.laplacian.ctimes \
	    --check-commute-scores-small test/$(1).smat.laplacian.commute-small \
	    --check-commute-scores-large test/$(1).smat.laplacian.commute-large
	$(2) ./graph_eigs test/$(1).smat -n -r -p > /dev/null
	python graph_eigs.py test/$(1).smat -n --check-eigs test/$(1).smat.normalized.eigs --check-ipar test/$(1).smat.normalized.ipar --check-resids test/$(1).smat.normalized.resids
	python graph_eigs.py test/$(1).smat --type=markov --check-eigs test/$(1).smat.normalized-markov.eigs --check-ipar test/$(1).smat.normalized-markov.ipar
	$(2) ./graph_eigs test/$(1).smat -m -r -p > /dev/null
	python graph_eigs.py test/$(1).smat -m --check-eigs test/$(1).smat.modularity.eigs --check-ipar test/$(1).smat.modularity.ipar --check-resids test/$(1).smat.modularity.resids
	
all_tests += test_$(1)
all_tests_clean += test/$(1).smat.*
endef
	

$(eval $(call TEST_template,tiny,))	
$(eval $(call TEST_template,karate,))	
$(eval $(call TEST_template,tapir,mpirun -np 4))
$(eval $(call TEST_template,polblogs-sym-cc,))
$(eval $(call TEST_template,Caltech36,mpirun -np 4))

test/element_iterator.o: scalapack_symmetric_eigen.cc
test/element_iterator: test/element_iterator.o 

test_element_iterator: test/element_iterator
	mpirun -np 4 test/element_iterator > /dev/null

all_tests_clean += test/element_iterator
all_tests += test_element_iterator

test: test_tiny test_karate test_element_iterator

all_tests: $(all_tests)

	
