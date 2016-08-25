
OPTFLAG = -O3
CXX = mpic++
FC = mpif77
CC = mpic++
FFLAGS += $(OPTFLAG) -funroll-loops
CXXFLAGS += -Wall -DBLACS_ALL -Iinclude -Wno-write-strings $(OPTFLAG)
LDLIBS += -Lscalapack/scalapack-2.0.2 -lscalapack -llapack -lblas -lgfortran

all : lapeigs graph_eigs identical_nodes

.PHONY : clean test clean_test

lapeigs.o :  scalapack_symmetric_eigen.hpp scalapack_symmetric_eigen.cc

lapeigs :  lapeigs.o mpiutil.o pdsyevr_tri.o

graph_eigs.o : graph_eigs_opts.hpp triplet_graph.hpp scalapack_symmetric_eigen.cc

graph_eigs : graph_eigs.o mpiutil.o pdsyevr_tri.o pdsyev_tri.o 

identical_nodes.o : triplet_graph.hpp

identical_nodes : identical_nodes.o 

clean:
	rm -rf graph_eigs graph_eigs.o lapeigs lapeigs.o mpiutil.o  pdsyevr_tri.o pdsyev_tri.o pdlawrite.o identical_nodes identical_nodes.o
	
clean_test:
	rm -rf $(all_tests_clean)

identical_nodes_test: identical_nodes
	rm -rf test/*.inodes
	./identical_nodes test/inodes_test1.smat
	diff test/inodes_test1.smat.inodes test/inodes_test1.correct
	./identical_nodes test/inodes_test2.smat
	diff test/inodes_test2.smat.inodes test/inodes_test2.correct
	./identical_nodes test/inodes_test3.smat test/inodes_test_three.smat.inodes
	diff test/inodes_test_three.smat.inodes test/inodes_test3.correct 

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
	    --check-commute-scores-large test/$(1).smat.laplacian.commute-large \
	    --check-pseudoinverse-diagonals test/$(1).smat.laplacian.pinvdiags \
	    --check-fiedler test/$(1).smat.laplacian.fiedler 
	$(2) ./graph_eigs test/$(1).smat -n -r -p > /dev/null
	python graph_eigs.py test/$(1).smat -n \
	    --check-eigs test/$(1).smat.normalized.eigs \
	    --check-ipar test/$(1).smat.normalized.ipar \
	    --check-resids test/$(1).smat.normalized.resids \
	    --check-fiedler test/$(1).smat.normalized.fiedler
	python graph_eigs.py test/$(1).smat --type=markov --check-eigs test/$(1).smat.normalized-markov.eigs --check-ipar test/$(1).smat.normalized-markov.ipar
	$(2) ./graph_eigs test/$(1).smat -m -r -p > /dev/null
	python graph_eigs.py test/$(1).smat -m --check-eigs test/$(1).smat.modularity.eigs --check-ipar test/$(1).smat.modularity.ipar --check-resids test/$(1).smat.modularity.resids
	
all_tests += test_$(1)
all_tests_clean += test/$(1).smat.*
endef
	

$(eval $(call TEST_template,tiny,))	
$(eval $(call TEST_template,karate,))	
$(eval $(call TEST_template,tapir,OMP_NUM_THREADS=1 mpirun -np 4))
$(eval $(call TEST_template,polblogs-sym-cc,))
#$(eval $(call TEST_template,Caltech36,mpirun -np 4))

test/element_iterator.o: scalapack_symmetric_eigen.cc
test/element_iterator: test/element_iterator.o 

test_element_iterator: test/element_iterator
	mpirun -np 4 test/element_iterator > /dev/null

all_tests_clean += test/element_iterator test/element_iterator.o test/*.inodes
all_tests += identical_nodes_test

test: test_tiny test_karate test_element_iterator identical_nodes_test

all_tests: $(all_tests)


