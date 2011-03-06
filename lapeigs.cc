/**
 * @file lapeigs.cc
 * Compute and save the eigenvalues of a Laplacian graph.
 * @author David F. Gleich
 * 
 * To compile
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <vector>

#include <mpi.h>
#include <blacs.h>
#include <scalapack.h>

#include "scalapack_symmetric_eigen.hpp"

#include "mpiutil.h"


const int lapeigs_version = 2;

void write_header(void) {
    mpi_world_printf("\n");
    for (int i=0; i<60; i++) {mpi_world_printf("%c", '-'); }
    mpi_world_printf("\n");
    mpi_world_printf("\n");
    mpi_world_printf("lapeigs\n");
    mpi_world_printf("=======\n");
    mpi_world_printf("\n");
    mpi_world_printf("Version %i\n", lapeigs_version);
    mpi_world_printf("David F. Gleich\n");
    mpi_world_printf("Copyright, 2010\n");
    mpi_world_printf("\n");
    for (int i=0; i<60; i++) {mpi_world_printf("%c", '-'); }
    mpi_world_printf("\n\n");
    
    int rank, size;
    char name[MPI_MAX_PROCESSOR_NAME+1];
    int namelen;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(name, &namelen);
    
    mpi_world_printf("nprocs = %i\n", size);
    
    mpi_world_sync_printf("[%3i] %s running...\n", rank, name);
}

bool write_vector(const char* filename, double *a, size_t n) {
    FILE *f = fopen(filename, "wt");
    if (f) {
        while (n > 0) {
            int val = fprintf(f,"%.18e\n",*a);
            if (val < 0) {
                // something failed
                // TODO make this output the processor information in COMM_WORLD
                printf("output failed\n");
                fclose(f);
                return (false);
            }
            n -= 1; 
            a += 1; 
        }
        fclose(f);
        return true;
    } else {
        printf("Could not open %s\n", filename);
        return false;
    }
}

/** Write a matrix in fortran format (column-oriented) */
bool write_matrix(const char* filename, double *a, size_t m, size_t n) {
    FILE *f = fopen(filename, "wt");
    if (f) {
        for (size_t i=0; i<m; ++i) {
            for (size_t j=0; j<n; ++j) {
                int val = fprintf(f,"%.18e ", a[i+j*m]);
                if (val < 0) {
                    // something failed
                    // TODO make this output the processor information in COMM_WORLD
                    printf("output failed\n");
                    fclose(f);
                    return (false);
                }
            }
            fprintf(f,"\n");
        }
        fclose(f);
        return true;
    } else {
        printf("Could not open %s\n", filename);
        return false;
    }
}

/** 
 * Usage:
 * 
 * main:triple_data g;
 * main:g.read_smat("myfile.smat");
 * main:g.distribute_data();
 */
struct triplet_data {
    int *r;
    int *c;
    double *v;
    int nrows;
    int ncols;
    int nnz;
    
    triplet_data() : r(NULL), c(NULL), v(NULL), nrows(0), ncols(0), nnz(0) {}
    ~triplet_data() {
        _free_smat();
    }
    
    void _free_smat() {
        if (r) { free(r); r = NULL; }
        if (c) { free(c); c = NULL; }
        if (v) { free(v); v = NULL; }
        nrows = 0;
        ncols = 0;
        nnz = 0;
    }
    
    void _alloc_smat(int _nr, int _nc, int _nz) {
        _free_smat();
        if ((_nr == 0 || _nc == 0) && _nz == 0) {
            // nothing to allocate
            return;
        }
        nrows = _nr;
        ncols = _nc;
        nnz = _nz;
        assert(nrows > 0);
        assert(ncols > 0);
        assert(nnz > 0);
        
        r = (int*)malloc(sizeof(int)*nnz);
        c = (int*)malloc(sizeof(int)*nnz);
        v = (double*)malloc(sizeof(double)*nnz);
        assert(r && c && v);
    }
        
    
    bool read_smat(const char* filename) {
        FILE *f = fopen(filename, "rt");
        if (f) {
            int m, n, nz;
            if (fscanf(f, "%i %i %i", &m, &n, &nz)==3) {
                _alloc_smat(m, n, nz);
                for (int nzi=0; nzi<nz; ++nzi) {
                    int i, j;
                    double a;
                    if (fscanf(f, "%i %i %lf", &i, &j, &a) != 3) {
                        fclose(f);
                        _free_smat();
                        return false;
                    }
                    if (i >= 0 && i < nrows && j>=0 && j<ncols) {
                        r[nzi] = i;
                        c[nzi] = j;
                        v[nzi] = a;
                    } else {
                        fclose(f);
                        _free_smat();
                        return false;
                    }
                }
            } else {
                fclose(f);
                return false;
            }
            fclose(f);
        } else {
            return false;
        }
        
        return true;
    }
    
    void broadcast(int ictxt) 
    {
        int nprow, npcol, myrow, mycol;
        Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
        int header[3];
        
        if (myrow==0 && mycol==0) {
            header[0] = nrows;
            header[1] = ncols;
            header[2] = nnz;
            
            Cigebs2d(ictxt, "All", "i-ring", 3, 1, header, 3);
            Cigebs2d(ictxt, "All", "i-ring", nnz, 1, r, nnz);
            Cigebs2d(ictxt, "All", "i-ring", nnz, 1, c, nnz);
            Cdgebs2d(ictxt, "All", "i-ring", nnz, 1, v, nnz);
        } else {
            Cigebr2d(ictxt, "All", "i-ring", 3, 1, header, 3, 0, 0);
            _alloc_smat(header[0], header[1], header[2]);
            Cigebr2d(ictxt, "All", "i-ring", nnz, 1, r, nnz, 0, 0);
            Cigebr2d(ictxt, "All", "i-ring", nnz, 1, c, nnz, 0, 0);
            Cdgebr2d(ictxt, "All", "i-ring", nnz, 1, v, nnz, 0, 0);
        }
    }
};

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

/**
 * The memory in A must be zeroed
 */
void assign_graph_modularity(int* desca, double* A, triplet_data& g)
{
    std::vector<int> degs(g.nrows);
    int izero = 0;
    int ictxt = desca[1], nrow=desca[2], ncol=desca[3],
      nbrow=desca[4], nbcol=desca[5];
    int myrow,mycol,nprow,npcol,ap, aq;
    
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
       
    
    assert(nrow==ncol);
    assert(nbrow==nbcol);
    
    ap = numroc_(&nrow, &nbrow, &myrow, &izero, &nprow);
    aq = numroc_(&ncol, &nbcol, &mycol, &izero, &npcol);
    
    // use fortran indexing here to avoid copy in indxl2g call
    for (int iloc=1; iloc <= ap; ++iloc) {
        for (int jloc=1; jloc <= aq; ++jloc) {
            int gi = indxl2g_(&iloc, &nbrow, &myrow, &izero, &nprow);
            int gj = indxl2g_(&jloc, &nbcol, &mycol, &izero, &npcol);
            A[(iloc-1) + (jloc-1)*ap] = 1./(gi + gj);
        }
    }
    
    

}

/** Initialize blacs
 * This is needed for many blacs implementations. 
 */
void blacs_init(void) {
    int iam, nprocs; // info for blacs init, should ignore.
    Cblacs_pinfo( &iam, &nprocs ) ;
}


/**
 * This is the main function on the BLACS computational grid.
 * There has been no BLACS init it.
 * 
 * This function assumes that all processors running it
 * are _in_ the blacs computational grid.
 */
int main_blacs(int argc, char **argv, int nprow, int npcol) {
    
    int ictxt;
    int myrow, mycol;
    triplet_data g;
    
    blacs_init();
    
    Cblacs_get(-1,0,&ictxt);
    Cblacs_gridinit(&ictxt, "R", npcol, npcol);    
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol );
    
    // check that this process is in the computational grid and
    // assert otherwise.
    assert(myrow < nprow);
    assert(mycol < npcol);
    
    if (myrow == 0 && mycol == 0) {
        if (argc < 3 || argc > 5) {
            printf("usage: lapeigs smatfile output [vectors=0|1|V] [minmemory=0|1]\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        const char* graphfilename = argv[1];
        printf("reading %s ...\n", graphfilename);
        if (g.read_smat(graphfilename) == false) {
            printf("Error reading %s.\n", graphfilename);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        printf("nrows = %i\n", g.nrows);
        printf("ncols = %i\n", g.ncols);
        printf("nnz = %i\n", g.nnz);
        
        printf("broadcasting ...\n");
    } 
    
    
    bool vectors = false;
    if (argc > 3) {
        if (strcmp(argv[3],"V")==0 || strcmp(argv[3],"1")==0) {
            vectors = true;
        }
    }
    bool minmemory = true;
    if (argc > 4) {
        if (strcmp(argv[4],"0")==0) {
            minmemory = false;
        }
    }
    
    g.broadcast(ictxt);
    
    printf("[%3i x %3i] received triplet data with %i non-zeros\n", 
        myrow, mycol, g.nnz);
    
    // wait here.
    Cblacs_barrier(ictxt, "A");
    
    // at this point, we have gotten all the data ready, so we 
    // enter the main eigenvalue routine
     
    // todo put this into a sub-routine
    {
        int n = g.nrows;
        // 176 seems like a good block size for a Nehalem, but
        // we need to check this number
        int nblock = std::min( 176, n/(2*nprow) );
        
        scalapack_distributed_matrix A;
        A.init(ictxt, n, n, nblock, nblock);
    
        // allocate local storage for the matrix
        printf("[%3i x %3i] allocating %Zu bytes for A (n=%i, nb=%i, ap=%i, aq=%i)\n",
            myrow, mycol, A.bytes(),
            A.n, A.nb, A.ap, A.aq );
            
        if (!A.allocate()) {
            printf("[%3i x %3i] couldn't allocate memory, exiting...", 
               myrow, mycol);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        
        assign_graph_laplacian(A.desc, A.A, g);
        
        Cblacs_barrier(ictxt, "A");
        
        scalapack_symmetric_eigen P(A);
        
        P.setup(vectors, minmemory, true, false);

        printf("[%3i x %3i] allocating %Zu bytes for eigenproblem\n",
            myrow, mycol, P.bytes());
        printf("[%3i x %3i] %i eigenvector elements; %i workspace\n",
            myrow, mycol, P.myZ*P.Z.ap*P.Z.aq, P.lwork);
        
        P.allocate();
    
        mpi_world_printf("Solving eigenproblem...\n");
        
        double t0;
        
        mpi_world_printf("Reducing to tridiagonal ...\n");
        t0 = MPI_Wtime();
        P.tridiag_reduce();
        double dt_red = MPI_Wtime() - t0;
            
        if (myrow==0 && mycol==0) {
            std::string str(argv[2]);
            str += ".tridiag";
            printf("Writing tridiagonal factors to %s...\n",
                str.c_str());
            if (write_matrix(str.c_str(), P.T, n, 2) == false) {
                printf("Writing failed... outputing to stdout\n");
                for (int i=0; i<n; ++i) {
                    printf("W(%i,1) = %.18e;\n", i+1, P.T[i]);
                    printf("W(%i,2) = %.18e;\n", i+1, P.T[i+n]);
                }   
            }
        }
        
        mpi_world_printf("Runing MRRR ...\n");
        t0 = MPI_Wtime();
        P.tridiag_compute();
        double dt_comp = MPI_Wtime() - t0;
        
        
        
        if (vectors) {
            std::vector<double> resids;
            P.residuals(resids);
            std::vector<double> ipars;
            P.Z.inverse_participation_ratios(ipars, true);
        }
        
        if (myrow==0 && mycol==0){
            printf("Computation:\n");
            printf(" vector = %i\n", (int)vectors);
            printf(" memory = %i\n", (int)minmemory);
            printf("      n = %d\n", n);
            printf("   grid = (%d,%d)\n",nprow,npcol);
            printf("  procs = %d procs\n", nprow*npcol);
            printf("  block = %d\n", nblock);
            printf(" rdtime = %.2fs\n", dt_red);
            printf(" cptime = %.2fs\n", dt_comp);
            printf("  ttime = %.2fs\n", dt_red+dt_comp);
        }
        
        if (myrow==0 && mycol==0) {
            if (write_vector(argv[2], P.values, (size_t)n) == false) {
                // error writint to argv[2], write to stdout
                printf("Writing eigenvalues to stdout\n");
                for (int i=0; i<n; ++i) {
                    printf("W(%i) = %.18e;\n", i+1, P.values[i]);
                }    
            }
        }
        
        P.free();
        A.free();
    }
    
    Cblacs_barrier(ictxt, "A");
    
    Cblacs_gridexit(ictxt);
    
    return 0;
}    
 
int main(int argc, char **argv) {
    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    write_header();
    
    // determine the number of processors to use in the square grid
    int np = (int)(floor(sqrt((double)size)));
    mpi_world_printf("using %i x %i grid of %i procs\n", np, np, np*np);
    
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
    if (rank < np*np) {
        main_blacs(argc, argv, np, np);
    } else {
        printf("[%3i] out of grid, exiting...\n", rank);
    }
    
    MPI_Finalize();
    return (0);
} 
