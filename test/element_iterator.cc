/**
 * @file element_iterator.cc
 * Testing routine for the element iterator
 */
 
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <mpi.h>
#include <blacs.h>
#include <scalapack.h>

#include "../scalapack_symmetric_eigen.hpp"
 
int main(int argc, char **argv) {
    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int np = (int)(floor(sqrt((double)size)));
    
    int ictxt;
    int myrow, mycol, nprow, npcol;
    
    int iam, nprocs; // info for blacs init, should ignore.
    Cblacs_pinfo( &iam, &nprocs ) ;
    
    Cblacs_get(-1,0,&ictxt);
    Cblacs_gridinit(&ictxt, "R", np, np);    
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol );
    
    bool root = (myrow == 0 && mycol == 0);
    
    scalapack_distributed_matrix A;
    int n = 25;
    int nb = 3;
    A.init(ictxt, n, n, nb, nb);
    bool rval = A.allocate();
    assert(rval);
        
    
    for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
            A.set(i, j, (double)(j+i*n));
        }
    }
    
    {
        if (root) {
            printf("creating iterator 1\n");
        }
        scalapack_distributed_matrix_element_iterator iter(A);
        assert(iter.root == root);
        if (root) {
            printf("calling next...\n");
        }
        while (iter.next()) {
            if (root) {
                printf("Got block (%i,%i) - (%i,%i)\n",
                    iter.istart, iter.jstart,
                    iter.istart+iter.isize - 1,
                    iter.jstart+iter.jsize - 1);
            }
            if (iter.root) {
                for (int j=0; j<iter.jsize; ++j) {
                    for (int i=0; i<iter.isize; ++i) {
                        printf("A(%i,%i) = %g\n",
                            iter.istart+i, iter.jstart+j, iter.element(i,j));
                    }
                }
            }
        }
    }
    
    {
        if (root) {
            printf("creating iterator 2\n");
        }
        scalapack_distributed_matrix_element_iterator iter(A, 1, 1, A.m);
        if (root) {
            printf("calling next...\n");
        }
        while (iter.next()) {
            if (root) {
                printf("Got block (%i,%i) - (%i,%i)\n",
                    iter.istart, iter.jstart,
                    iter.istart+iter.isize - 1,
                    iter.jstart+iter.jsize - 1);
            }
            if (iter.root) {
                for (int j=0; j<iter.jsize; ++j) {
                    for (int i=0; i<iter.isize; ++i) {
                        printf("A(%i,%i) = %g\n",
                            iter.istart+i, iter.jstart+j, iter.element(i,j));
                    }
                }
            }
        }
    }
    
    {
        if (root) {
            printf("creating iterator 3\n");
        }
        scalapack_distributed_matrix_element_iterator iter(A, 1, 1, 2*A.m);
        if (root) {
            printf("calling next...\n");
        }
        while (iter.next()) {
            if (root) {
                printf("Got block (%i,%i) - (%i,%i)\n",
                    iter.istart, iter.jstart,
                    iter.istart+iter.isize - 1,
                    iter.jstart+iter.jsize - 1);
            }
            if (iter.root) {
                for (int j=0; j<iter.jsize; ++j) {
                    for (int i=0; i<iter.isize; ++i) {
                        printf("A(%i,%i) = %g\n",
                            iter.istart+i, iter.jstart+j, iter.element(i,j));
                    }
                }
            }
        }
    }
    
    {
        if (root) {
            printf("creating iterator 4\n");
        }
        scalapack_distributed_matrix_element_iterator iter(A, 1, 1, A.n*A.m);
        if (root) {
            printf("calling next...\n");
        }
        while (iter.next()) {
            if (root) {
                printf("Got block (%i,%i) - (%i,%i)\n",
                    iter.istart, iter.jstart,
                    iter.istart+iter.isize - 1,
                    iter.jstart+iter.jsize - 1);
            }
            if (iter.root && root) {
                for (int j=0; j<iter.jsize; ++j) {
                    for (int i=0; i<iter.isize; ++i) {
                        printf("A(%i,%i) = %g\n",
                            iter.istart+i, iter.jstart+j, iter.element(i,j));
                    }
                }
            }
        }
    }
    
    MPI_Finalize();
    
    return 0;
}
