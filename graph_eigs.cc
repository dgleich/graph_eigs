/**
 * @file graph_eigs.cc
 * Compute all the eigenvalues (and possibly eigenvectors) of a large graph.
 */

/** 
 * David F. Gleich
 * Copyright, 2010
 */

/** 
 * History
 * -------
 * :2010-10-07: Initial coding
 * :2010-11-10: resuming coding after finishing lapeigs in C with all
 *   desired scalapack features.
 * 
 * Todo
 * ----
 * TODO fix adjacency/laplacian/nlaplacian to handle repeated edges
 * TODO add report with output:
 *   - graph name
 *   - number of vertices
 *   - number of edges
 *   - matrix type
 *   - ... (more options)
 *   - nprocs
 *   - grid size
 *   - output :
 *      - eigenvectors
 *      - tridiag
 *      ...
 *   - timing:
 *      - graph load
 *      - assign
 *      - tridiag
 *      - eigensolve
 *      - residuals
 *      - iparscores
 *      - total 
 *   - memory:
 *      - matrix
 *      - eigensystem
 *      - outputs
 *   - eigenvalues:
 *      - largest
 *      - smallest
 *      - median
 *      - mode
 *   - residuals
 *      - maximum
 *      - large
 *      (other useful info on the residuals)
 *   - iparscores
 *      - maximum
 *      (other useful info 
 * 
 * The entire report shouldn't be too large on disk, and should
 * be in a human and computer readerable format (json, yaml)
 */

#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>

#include <vector>

#include <mpi.h>
#include <blacs.h>
#include <scalapack.h>

#include "scalapack_symmetric_eigen.hpp"

#include "mpiutil.h"

#include "triplet_graph.hpp"
#include "graph_eigs_opts.hpp"

const int graph_eigs_version = 2;

void write_header(void) {
    mpi_world_printf("\n");
    for (int i=0; i<60; i++) {mpi_world_printf("%c", '-'); }
    mpi_world_printf("\n");
    mpi_world_printf("\n");
    mpi_world_printf("graph_eigs\n");
    mpi_world_printf("==========\n");
    mpi_world_printf("\n");
    mpi_world_printf("Version %i\n", graph_eigs_version);
    mpi_world_printf("David F. Gleich\n");
    mpi_world_printf("Copyright, 2010\n");
    mpi_world_printf("\n");
    for (int i=0; i<60; i++) {mpi_world_printf("%c", '-'); }
    mpi_world_printf("\n\n");
    
    if (opts.verbose) {
        int rank, size;
        char name[MPI_MAX_PROCESSOR_NAME+1];
        int namelen;
        
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Get_processor_name(name, &namelen);
        
        mpi_world_printf("nprocs = %i\n", size);

        mpi_world_sync_printf("[%3i] %s running...\n", rank, name);
    }
}

void broadcast_triplet_data(int ictxt, triplet_data& graph ) 
{
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
    int header[3];
    
    if (myrow==0 && mycol==0) {
        header[0] = graph.nrows;
        header[1] = graph.ncols;
        header[2] = graph.nnz;
        
        Cigebs2d(ictxt, "All", "i-ring", 3, 1, header, 3);
        Cigebs2d(ictxt, "All", "i-ring", graph.nnz, 1, graph.r, graph.nnz);
        Cigebs2d(ictxt, "All", "i-ring", graph.nnz, 1, graph.c, graph.nnz);
        Cdgebs2d(ictxt, "All", "i-ring", graph.nnz, 1, graph.v, graph.nnz);
    } else {
        Cigebr2d(ictxt, "All", "i-ring", 3, 1, header, 3, 0, 0);
        graph._alloc_data(header[0], header[1], header[2]);
        Cigebr2d(ictxt, "All", "i-ring", graph.nnz, 1, graph.r, graph.nnz, 0, 0);
        Cigebr2d(ictxt, "All", "i-ring", graph.nnz, 1, graph.c, graph.nnz, 0, 0);
        Cdgebr2d(ictxt, "All", "i-ring", graph.nnz, 1, graph.v, graph.nnz, 0, 0);
        if (opts.verbose) {
            printf("[%3i x %3i] received triplet data with %i edges\n", 
                myrow, mycol, graph.nnz);
        }
    }
}

/** Assignment routines for different matrices.
 */

/**
 * The memory in A must be zeroed
 */
void assign_graph_adjacency(int* desca, double* A, triplet_data& g)
{
    std::vector<int> degs(g.nrows);
    
    for (int nzi=0; nzi<g.nnz; ++nzi) {
        int i = g.r[nzi]+1;
        int j = g.c[nzi]+1;
        double a = 1.;
        pdelset_(A, &i, &j, desca, &a);
    }
}

/**
 * The memory in A must be zeroed
 */
void assign_graph_laplacian(int* desca, double* A, triplet_data& g)
{
    std::vector<int> degs(g.nrows);
    
    // here, we use the function pdelset_
    // which only performs an assignment on the correct node
    // so we don't have to worry about the block-cyclic organization
    
    for (int nzi=0; nzi<g.nnz; ++nzi) {
        int i = g.r[nzi]+1;
        int j = g.c[nzi]+1;
        double a = -1.;
        if (i != j) {
            degs[g.r[nzi]] += 1;
            pdelset_(A, &i, &j, desca, &a);
        }
    }
    
    for (int i=1; i<=g.nrows; ++i) {
        double a = degs[i-1];
        pdelset_(A, &i, &i, desca, &a);
    }

}

void assign_graph_normalized_laplacian(int* desca, double* A, triplet_data& g)
{
    std::vector<int> degs(g.nrows);
    std::vector<int> diag(g.nrows,0.);
    
    // here, we use the function pdelset_
    // which only performs an assignment on the correct node
    // so we don't have to worry about the block-cyclic organization
    
    for (int nzi=0; nzi<g.nnz; ++nzi) {
        degs[g.r[nzi]] += 1;
        int i = g.r[nzi]+1;
        int j = g.c[nzi]+1;
        if (i == j) {
            diag[i-1] += 1;
        }
    }
    
    for (int nzi=0; nzi<g.nnz; ++nzi) {
        int i = g.r[nzi]+1;
        int j = g.c[nzi]+1;
        double a = -1.;
        if (degs[i-1] > 0 && degs[j-1] > 0) {
            a = -1./(sqrt((double)degs[i-1])*sqrt((double)degs[j-1]));
        } else {
            // we'll have to figure out how to handle this case later
            assert(false);
        }
        if (i != j) {
            pdelset_(A, &i, &j, desca, &a);
        } 
    }
    
    for (int i=1; i<=g.nrows; ++i) {
        double a = 1.0 - (double)diag[i-1]/(double)degs[i-1];
        pdelset_(A, &i, &i, desca, &a);
    }

}

/** Initialize blacs
 * This is needed for many blacs implementations. 
 */
void blacs_init(void) {
    int iam, nprocs; // info for blacs init, should ignore.
    Cblacs_pinfo( &iam, &nprocs ) ;
}


int main_blacs(int argc, char **argv, int nprow, int npcol) 
{
    int ictxt;
    int myrow, mycol;
    
    blacs_init();
    
    Cblacs_get(-1,0,&ictxt);
    Cblacs_gridinit(&ictxt, "R", npcol, npcol);    
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol );
    
    bool root = (myrow == 0 && mycol == 0);
    
    // check that this process is in the computational grid and assert otherwise.
    assert(myrow < nprow);
    assert(mycol < npcol);
    
    triplet_data g;
    
    if (root) {
        const char* gfile = opts.graph_filename.c_str();
        printf("Loading graph %s.\n", gfile);
        if (g.read_smat(gfile) == false) {
            printf("Error reading %s.\n", gfile);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        printf("\n");
        printf("Graph information\n");
        printf("  Vertices: %i\n", g.nrows);
        printf("  Edges: %i\n", g.nnz);
        printf("\n");
        
        if (opts.verbose) {
            printf("Broadcasting graph data.\n");
        }
        
    }
    
    broadcast_triplet_data(ictxt, g);
    
    // we now allocate memory for the dense matrix to store the graph
    scalapack_distributed_matrix A;
    int n = g.nrows;
    A.init(ictxt, n, n, opts.nb, opts.nb);
    
    // allocate local storage for the matrix
    if (opts.verbose) {
        printf("[%3i x %3i] allocating %Zu bytes for A (n=%i, nb=%i, ap=%i, aq=%i)\n",
            myrow, mycol, A.bytes(),
            A.n, A.nb, A.ap, A.aq );
    }
        
    if (!A.allocate()) {
        printf("[%3i x %3i] couldn't allocate memory, exiting...", 
           myrow, mycol);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    switch (opts.matrix) {
        case graph_eigs_options::adjacency_matrix:
            //assign_graph_laplacian(A.desc, A.A, g);
            mpi_world_printf("Constructing the adjacency matrix.\n");
            assign_graph_adjacency(A.desc, A.A, g);
            break;
            
        case graph_eigs_options::laplacian_matrix:
            mpi_world_printf("Constructing the Laplacian matrix.\n");
            assign_graph_laplacian(A.desc, A.A, g);
            break;
            
        case graph_eigs_options::normalized_laplacian_matrix:
            mpi_world_printf("Constructing the normalized Laplacian matrix.\n");
            assign_graph_normalized_laplacian(A.desc, A.A, g);
            break;
            
        case graph_eigs_options::modularity_matrix:
            mpi_world_printf("Constructing the modularity matrix.\n");
            assert(false);
            break;
    }
    
    scalapack_symmetric_eigen P(A);
        
    P.setup(opts.eigenvectors, opts.minmemory, opts.tridiag);
    
    return (0);
}


int main(int argc, char **argv) {
    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int tocontinue;
    if (rank == 0) {
        tocontinue = (int)parse_command_line_arguments(argc, argv);
        if (tocontinue) {
            opts.setup(); // setup things like filenames
            tocontinue = opts.check_filenames();
        }
        
        MPI_Bcast(&tocontinue, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        MPI_Bcast(&tocontinue, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
       
    if (!tocontinue) {
        MPI_Finalize();
        return (-1);
    }
    
    opts.distribute(); 
    
    write_header();
    
    
    
    // determine the number of processors to use in the square grid
    int np = (int)(floor(sqrt((double)size)));
    if (opts.verbose) {
        mpi_world_printf("using %i x %i grid of %i procs\n", np, np, np*np);
    }
    
    #ifdef BLACS_ALL
    if (size != np*np) {
        mpi_world_printf("BLACS_ALL defined but BLACS grid=%i and MPI size=%i\n",
            np*np, size);
        mpi_world_printf("BLACS_ALL requires all MPI procs to be in the BLACS grid\n");
        mpi_world_printf("exiting...\n");
        MPI_Finalize();
        return (-1);
    }
    #endif /* BLACS_ALL */
    

    // initialize the blacs grid
    int rval = 0;
    if (rank < np*np) {
        rval = main_blacs(argc, argv, np, np);
    } else {
        if (opts.verbose) {
            printf("[%3i] out of grid, exiting...\n", rank);
        }
    }
    
    MPI_Finalize();
    
    return (rval);
}
