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
 * :2010-11-14: Added Markov matrix
 * :2010-11-17: Added graph checking
 * :2011-02-15: Added commute time computation
 * :2011-03-05: Added commute time score output
 * :2011-03-06: Added Fiedler computation
 * :2011-03-08: Added pseudo-inverse diagonal output
 * 
 * 
 * Todo
 * ----
 * TODO output extremal eigenvalues for checking
 * TODO add options for pseudoinverse of the diagonal
 * TODO Add switching eigensolver
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
 * The entire report shouldn't be too large on disk, and should
 * be in a human and computer readerable format (json, yaml)
 */

#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>

#include <string>
#include <vector>

#include <mpi.h>
#include <blacs.h>
#include <scalapack.h>

#include "scalapack_symmetric_eigen.hpp"

#include "mpiutil.h"

#include "triplet_graph.hpp"
#include "graph_eigs_opts.hpp"
#include "permutation.hpp"

const int graph_eigs_version = 4;

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

struct wall_time_list {
    typedef std::pair< std::string, double > event;
    std::vector< event > events;
    
    double t0;
    double ttime;
    
    void start_event( const char* name ) {
        event e;
        e.first = std::string(name);
        e.second = -1;
        events.push_back(e);
        t0 = MPI_Wtime();
    }
    
    void end_event() {
        assert(events.size() > 0);
        events.back().second = MPI_Wtime() - t0;
        ttime += events.back().second;
    }
    
    void report_event( unsigned int indent = 0 ) {
        assert(events.size() > 0);
        event last = events.back();
        while (indent > 0) {
            printf(" ");
            indent -= 1;
        }
        
        double rtime = last.second;
        if (rtime < 0) {
            rtime = MPI_Wtime() - t0;
        }
        
        printf("%s time: %.1lf sec.\n", last.first.c_str(), rtime); 
    }
};

wall_time_list tlist;

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
        printf("Error: Could not open %s\n", filename);
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
        printf("Error: Could not open %s\n", filename);
        return false;
    }
}

/** Assignment routines for different matrices.
 */

/**
 * The memory in A must be zeroed
 */
void assign_graph_adjacency(scalapack_distributed_matrix& A, triplet_data& g)
{
    A.set_to_constant(0.);
    
    for (int nzi=0; nzi<g.nnz; ++nzi) {
        A.incr(g.r[nzi], g.c[nzi], 1.);
    }
}

/**
 */
void assign_graph_laplacian(scalapack_distributed_matrix& A, triplet_data& g)
{
    std::vector<int> degs(g.nrows, 0);
    
    A.set_to_constant(0.);
    
    // the diagonal entries of the laplacian cancel out, so don't count them.
    
    for (int nzi=0; nzi<g.nnz; ++nzi) {
        int i = g.r[nzi];
        int j = g.c[nzi];
        A.incr(i,j,-1.);
        degs[g.r[nzi]] += 1;
    }
    
    for (int i=0; i<g.nrows; ++i) {
        A.incr(i,i,(double)degs[i]);
    }
}

/**
 * This routine creates a non-symmetric matrix.  It cannot be used
 * to construct a matrix for scalapack_symmetric_eigen.  Instead, 
 * see how we use the normalized laplacian to construct the Markov matrix.
 */
void assign_graph_markov(scalapack_distributed_matrix& A, triplet_data& g)
{
    std::vector<int> degs(g.nrows, 0);
    
    A.set_to_constant(0.);
    
    // the diagonal entries of the laplacian cancel out, so don't count them.
    for (int nzi=0; nzi<g.nnz; ++nzi) {
        degs[g.r[nzi]] += 1;
    }
    
    for (int nzi=0; nzi<g.nnz; ++nzi) {
        int i = g.r[nzi];
        int j = g.c[nzi];
        assert(degs[i] > 0);
        assert(degs[j] > 0);
        double a = 1./((double)degs[i]);
        
        A.incr(g.r[nzi],g.c[nzi],a);
    }
  
}

void assign_graph_normalized_laplacian(scalapack_distributed_matrix& A, triplet_data& g)
{
    std::vector<int> degs(g.nrows,0.);
    
    A.set_to_constant(0.);
    
    for (int nzi=0; nzi<g.nnz; ++nzi) {
        degs[g.r[nzi]] += 1;
    }
    
    for (int nzi=0; nzi<g.nnz; ++nzi) {
        int i = g.r[nzi];
        int j = g.c[nzi];
        assert(degs[i] > 0);
        assert(degs[j] > 0);
        double a = -1./(sqrt((double)degs[i])*sqrt((double)degs[j]));
        
        A.incr(g.r[nzi],g.c[nzi],a);
    }
    for (int i=0; i<g.nrows; i++) {
        A.incr(i,i,1.0);
    }
}

void assign_graph_modularity(scalapack_distributed_matrix& A, triplet_data& g)
{
    std::vector<int> degs(g.nrows,0.);
    
    A.set_to_constant(0.);
    
    double vol = 0;
    for (int nzi=0; nzi<g.nnz; ++nzi) {
        degs[g.r[nzi]] += 1;
        vol += 1.;
    }
    
    double ivol = 1./vol;
    
    for (int li = 0; li<A.ap; ++li) {
        for (int lj = 0; lj<A.aq; ++lj) {
            int gi, gj;
            A.local2global(li,lj,gi,gj);
            A.A[li+lj*A.ap] = -1.0*ivol*(double)degs[gi]*(double)degs[gj];
        }
    }
    
    for (int nzi=0; nzi<g.nnz; ++nzi) {
        A.incr(g.r[nzi], g.c[nzi], 1.);
    }
    
}

void assign_graph_matrix(triplet_data& g, scalapack_distributed_matrix& A, 
    graph_eigs_options::matrix_type t)
{
    switch (t) {
        case graph_eigs_options::adjacency_matrix:
            mpi_world_printf("Constructing the adjacency matrix.\n");
            assign_graph_adjacency(A, g);
            break;
            
        case graph_eigs_options::laplacian_matrix:
            mpi_world_printf("Constructing the Laplacian matrix.\n");
            assign_graph_laplacian(A, g);
            break;
            
        case graph_eigs_options::normalized_laplacian_matrix:
            mpi_world_printf("Constructing the normalized Laplacian matrix.\n");
            assign_graph_normalized_laplacian(A, g);
            break;
            
        case graph_eigs_options::modularity_matrix:
            mpi_world_printf("Constructing the modularity matrix.\n");
            assign_graph_modularity(A, g);
            break;
    }
}

/** Initialize blacs
 * This is needed for many blacs implementations. 
 */
void blacs_init(void) {
    int iam, nprocs; // info for blacs init, should ignore.
    Cblacs_pinfo( &iam, &nprocs ) ;
}

/** This function tries to write to a file, but then hits stdout otherwise.
 * 
 * if matrix is false, then n is ignored.
 */
void write_data_safely(const char *type, const char* vname, unsigned int indent, 
        const char *filename, double* a, size_t m, size_t n, bool matrix) 
{
    while (indent > 0) {
        printf(" ");
        indent -= 1;
    }
    printf("  Writing %s to %s.\n", type, filename);
    if (matrix) {
        if (write_matrix(filename, a, m, n) == false) {
            printf("Error: writing %s to stdout\n", type);
            for (size_t i=0; i<m; ++i) {
                for (size_t j=0; j<n; ++j) {
                    printf("%s(%Zu,%Zu) = %.18e;\n", vname, i+1, j+1, a[i + j*m]);
                }
            }    
        }
    } else {
        if (write_vector(filename, a, (size_t)m) == false) {
            printf("Error: writing %s to stdout\n", type);
            for (size_t i=0; i<m; ++i) {
                printf("%s(%Zu) = %.18e;\n", vname, i+1, a[i]);
            }    
        }
    }
}

/** Output any residuals that are larger than a set scale. 
 * This function is useful to indicate any problems with the computation.
 * @param resids the residuals associated with each eigenvector.
 * @param scale a 
 * @param tol the large residual tolerance
 */
void find_large_residuals(std::vector<double>& resids, double *scale, double tol)
{
    printf("  Checking for residuals larger than %.2e.\n", tol);
    for (size_t i=0; i<resids.size(); ++i) {
        double err = fabs(resids[i])/(fabs(scale[i])+1.);
        if (err > tol) {
            if (err > 1000*tol) {
                printf(" **** Residual %Zu is EXTREMELY large: %.18e ; eval=%.2e\n", 
                    i, err, scale[i]);
            } else {
                printf("      Residual %Zu is large: %.18e ; eval=%.2e\n",
                    i, err, scale[i]);
            }
        }
    }
}

FILE* fopen_with_error(const char* filename, const char* mode) {
    FILE *f = fopen(filename, mode);
    if (!f) {
        fprintf(stderr, "Error: Could not open %s with mode %s\n", 
            filename, mode);
    }
    return f;
}

/** Output the largest and smallest elements in a matrix to files.
 * 
 * Given a matrix, for each column, we output the `nelem` largest
 * and smallest elements to `large_filename` and `small_filename`,
 * respectively.  This involves saving each column on one processor
 * and then sorting the elements, and outputing the respective
 * sets to each file.  The file format is:
 * File: Header\n Element*[<nelem>*<ncols>]\n
 * Header: <nrows> <ncols> <nelem>*<ncols>
 * Element: <rowid> <colid> <value>\n
 */
void write_largest_and_smallest_matrix_column_elements(
    scalapack_distributed_matrix& A,
    int nelem,
    std::string large_filename, std::string small_filename)
{
    // check that nelem isn't too large
    if (nelem > A.m) { nelem = A.m; }
    
    // create an iterator which gets a column of the matrix at once,
    // with a A.m sized buffer
    scalapack_distributed_matrix_element_iterator iter(A, 1, 1, A.m);
    
    //printf("[%3i x %3i] Starting iterator...\n",
    //     A.myrow, A.mycol);
    
    FILE *sfile=NULL, *lfile=NULL;
    if (iter.root) {
        sfile = fopen_with_error(small_filename.c_str(), "wt");
        lfile = fopen_with_error(large_filename.c_str(), "wt");
        if (!sfile || !lfile) {
            return;
        }
        
        // write the headers
        fprintf(sfile, "%i %i %i\n", A.m, A.n, nelem*A.n);
        fprintf(lfile, "%i %i %i\n", A.m, A.n, nelem*A.n);
    }    
    
    while (iter.next()) {
        
        // this will get one column
        if (iter.root) {
            //printf("[%3i x %3i] Next column (%i,%i) - (%i,%i)\n",
            //    A.myrow, A.mycol, iter.istart, iter.jstart, iter.isize, iter.jsize);
            assert(iter.isize == A.m);
            assert(iter.jsize == 1);
            std::vector<int> order(0);
            sort_permutation(iter.work, order);
            // the array is now sorted in decreasing order
            // the largest elements are in iter.work[0:nelem]
            // the smallest element are in iter.work[end-nelem+1:]
            
            assert(nelem <= (int)order.size());
            
            // output the largest
            for (int i=0; i<nelem; ++i) {
                fprintf(lfile, "%i %i %.18g\n", 
                    order[i], iter.jstart-1, iter.work[order[i]]);
            }
            // output the smallest elements, except for the last element,
            // which is always the diagonal one.
            for (int i=0; i<nelem; ++i) {
                int ei = order.size() - nelem + i;
                fprintf(sfile, "%i %i %.18g\n",
                    order[ei], iter.jstart-1, iter.work[order[ei]]);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    if (iter.root) {
        fclose(sfile);
        fclose(lfile);
    }
}    

/** Compute and output data on commute times.
 * 
 * This uses extra work memory used for computing the eigenvalues.
 * 
 */
template <typename EigenSolver>
void compute_commute_times(scalapack_distributed_matrix& A, EigenSolver &P) {
    bool root = A.myrow == 0 && A.mycol == 0;
    if (root) {
        printf("  Computing commute-time scores\n");
    }
    scalapack_distributed_matrix C;
    if (root) {
        tlist.start_event("commute");
    }
    P.pseudoinverse(C); // compute the pseudo-inverse
    std::vector<double> diags(A.n,0.);
    C.diagonals(&diags[0]);
    
    int gi, gj;
    for (int j=0; j<C.aq; ++j) {
        for (int i=0; i<C.ap; ++i) {
            C.local2global(i, j, gi, gj);
            // compute the commute times
            C.A[i+j*C.ap] = diags[gi] + diags[gj] - 2*C.A[i+j*C.ap];
            if (gi==gj) {
                C.A[i+j*C.ap] = 0.;
            }
        }
    }
    if (root) {
        tlist.report_event(6);
    }
    
    if (opts.pseudoinverse_diagonals) {
        if (root) {
            write_data_safely("pseudoinverse diagonal elements", "pid", 2,
                opts.pseudoinverse_diagonals_filename.c_str(), 
                &diags[0], A.n, 1, false);
        }
    }
    
    if (opts.commute_scores) {
        if (root) {
            printf("    Writing %i largest commute-time elements to %s.\n",
                opts.ncommute_scores, opts.commute_large_scores_filename.c_str());
            printf("    Writing %i smallest commute-time elements to %s.\n",
                opts.ncommute_scores, opts.commute_small_scores_filename.c_str());
        }
        write_largest_and_smallest_matrix_column_elements(C, 
            opts.ncommute_scores,
            opts.commute_large_scores_filename, 
            opts.commute_small_scores_filename);
    }
    
    if (opts.commute_all) {
        if (root) {
            printf("    Writing commute-time matrix to %s.\n",
                    opts.commute_all_filename.c_str());
        }
        C.write(opts.commute_all_filename.c_str(), 0, 0);
    }
}    

/** Find the fielder vector
 * The fiedler vector is the eigenvector for the first non-zero eigenvalue
 * of a Laplacian matrix.  
 * 
 * This function works with the absolute values of the eigenvalues to 
 * avoid issues with slightly negative eigenvalues exceeding the 
 * a tolerance.
 * 
 * An eigenvalue is considered zero if fabs(lambda) < n*sigmamax*2.2e-16
 * 
 * @param fielder 
 * @param index an output variable for the index of the Fiedler vector
 * @return true if the Fiedler vector was found.
 */
bool find_fiedler(double* values, scalapack_distributed_matrix& Z,
    double *fiedler, int* index)
{  
    // this function does not assume that the eigenvalues are sorted.
      
    // first determine the zero eigenvalue tolerance
    // find the largest eigenvalue value of A
    double sigmamax = values[0];
    int imax = 0;
    for (int i=1; i<Z.n; ++i) {
        if (fabs(values[i]) > sigmamax) {
            sigmamax = fabs(values[i]);
            imax = i;
        } 
    }
    // set the zero tolerance
    double tol = (double)Z.n*2.2e-16*sigmamax;
    
    // now find the smallest non-zero eigenvalue
    if (sigmamax < 2.2e-16*(double)Z.n) {
        // there are no non-zero eigenvalues...
        return false;
    }
    
    double nzmin = sigmamax;
    int nzind = imax;
    for (int i=0; i<Z.n; ++i) {
        if (values[i] > tol && values[i] < nzmin) {
            nzmin = values[i];
            nzind = i;
        }
    }
    // get the column
    if (fiedler) {
        Z.column(nzind, fiedler);
    }
    
    if (index) {
        *index = nzind;
    }
    
    return true;
}




/** Compute data on the Markov matrix from the normalized Laplacian.
 * For an undirected graph, the eigenvalues and vectors of the 
 * Markov chain transition matrix can be computed via the eigenvalues
 * and eigenvectors of the normalized Laplacian matrix.  This function
 * transforms the eigenvectors and values, and also outputs the
 * new information.
 *
 * *IT DESTROYS INFORMATION IN A and P*, the matrix and problem.
 * 
 * @param g the graph
 * @param A the current matrix
 * @param P the current eigensystem
 * */
void output_markov_data(triplet_data& g,
    scalapack_distributed_matrix& A, 
    scalapack_symmetric_eigen& P)
{
    int myrow, mycol, nprow, npcol;
    Cblacs_gridinfo(A.ictxt, &nprow, &npcol, &myrow, &mycol );
    bool root = myrow == 0 && mycol == 0;
    
    assert(opts.matrix == opts.normalized_laplacian_matrix);
        
    // change the eigenvalues
    for (int i=0; i<P.n; ++i) {
        P.values[i] = 1.0 - P.values[i];
    }
    
    if (opts.eigenvalues && root)  {
        write_data_safely("Markov eigenvalues", "mw", 2,
            opts.markov_values_filename.c_str(), P.values, P.n, 1, false);
    }
        
    if (opts.eigenvectors) {
        // compute graph degrees
        assert(g.nrows == P.n);
        std::vector<double> degs(g.nrows, 0.);
        for (int nzi=0; nzi<g.nnz; ++nzi) {
            degs[g.r[nzi]] += 1.;
        }
        // invert graph degres
        for (int i=0; i<g.nrows; ++i) {
            assert(degs[i] > 0);
            degs[i] = sqrt(1.0/(double)degs[i]);
        }
        
        // scale the eigenvectors
        P.Z.scale_rows(&degs[0], 1.0);
        
        if (opts.residuals) {
            std::vector<double> resids;
            tlist.start_event("markov_residuals");
            assign_graph_markov(P.A, g);
            P.reloaded_residuals(resids);
            tlist.end_event();
            if (root) { tlist.report_event(6); }
            
            if (root) {
                write_data_safely("Markov residuals", "mr", 2,
                    opts.markov_residuals_filename.c_str(), &resids[0], P.n, 1, false);
            }
            
            if (root) {
                find_large_residuals(resids, P.values, 2.e-16*P.n);
            }
        }
        if (opts.iparscores) {
            std::vector<double> ipars;
            tlist.start_event("markov_iparscores");
            P.Z.inverse_participation_ratios(ipars, false);  
            tlist.end_event();
            if (root) { tlist.report_event(6); }
            
            if (root) {
                write_data_safely("Markov iparscores", "mp", 2,
                    opts.markov_ipar_filename.c_str(), &ipars[0], P.n, 1, false);
            }
        }
        if (opts.vectors) {
            // write the matrix on the root processor
            if (root) {
                printf("    Writing eigenvectors to %s.\n", 
                    opts.markov_vectors_filename.c_str());
            }
            P.Z.write(opts.markov_vectors_filename, 0, 0);
        }
    }
}   

template <typename EigenSolver>
int main_compute(bool root, triplet_data& g, scalapack_distributed_matrix& A) 
{
    int nprow, npcol, myrow, mycol;
    int n = A.n;
    Cblacs_gridinfo(A.ictxt, &nprow, &npcol, &myrow, &mycol);
    
    EigenSolver P(A);
    
    bool extra = false;
    if (opts.matrix == opts.laplacian_matrix &&
        opts.commute && (opts.commute_all || opts.commute_scores)) {
        extra = true;
    }
    
    P.setup(opts.eigenvectors, opts.minmemory, opts.tridiag, extra);
    
    if (opts.verbose) {
        printf("[%3i x %3i] allocating %Zu bytes for eigenproblem; "
               "Z=%i, work=%i, iwork=%i\n",
            myrow, mycol, P.bytes(), P.vectors*P.Z.ap*P.Z.aq, P.lwork, 
            P.liwork);
    }
    
    if (!P.allocate()) {
        printf("[%3i x %3i] couldn't allocate memory, exiting...", 
           myrow, mycol);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    if (opts.eigenvectors == false && opts.eigenvalues == false) {
        mpi_world_printf("Computing tridiagonal redunction only.\n");
        tlist.start_event("tridiag_reduce");
        P.tridiag_reduce();
        tlist.end_event();
        if (root) { tlist.report_event(2); }
    } else {
        if (opts.eigenvectors == false) {
            mpi_world_printf("Computing eigenvalues.\n");
        } else {
            mpi_world_printf("Computing eigenvalues and eigenvectors.\n");
        }

        // this is for debugging issues with pdsyevr_tri
        // mpi_world_printf("  Full eigensolve.\n");
        //tlist.start_event("full eigensolve");
        //P.compute();
        //tlist.end_event();
        
        
        
            
        mpi_world_printf("  Reducing to tridiagonal.\n");
        
        tlist.start_event("tridiag_reduce");
        P.tridiag_reduce();
        tlist.end_event();
        if (root) { tlist.report_event(4); }
        
        if (opts.tridiag && myrow==0 && mycol==0)  {
            write_data_safely("tridiagonal factors", "T", 2,
                opts.tridiag_filename.c_str(), P.T, n, 2, true);
        }
        
        if (opts.eigenvectors == false) {
            mpi_world_printf("  Finding eigenvalues.\n");    
        } else {
            mpi_world_printf("  Finding eigenvalues and eigenvectors.\n");
        }
        
        tlist.start_event("eigensolve");
        P.tridiag_compute();
        tlist.end_event();
        if (root) { tlist.report_event(4); } 
        
        if (opts.eigenvalues && root)  {
            write_data_safely("eigenvalues", "W", 2,
                opts.values_filename.c_str(), P.values, n, 1, false);
        }
        
        if (opts.eigenvectors) {
            
            if (opts.vectors) {
                // write the matrix on the root processor
                mpi_world_printf("    Writing eigenvectors to %s.\n",
                    opts.vectors_filename.c_str());
                P.Z.write(opts.vectors_filename, 0, 0);
                //write_matrix("test.matrix", P.Z.A, P.Z.aq, P.Z.ap);
            }
            
            // output the Fiedler vector
            // TODO make this a function
            if ((opts.matrix == opts.normalized_laplacian_matrix ||
                opts.matrix == opts.laplacian_matrix) && opts.fiedler) {
                std::vector<double> fiedler(P.n, 0.);
                int fvind = 0;
                bool rval = find_fiedler(P.values, P.Z, &fiedler[0], &fvind);
                
                if (rval == false) {
                    mpi_world_printf("  Could not find fiedler vector.\n");
                } else {
                    double lam = P.values[fvind];
                    double gap0 = 1./0.;
                    double gap1 = 1./0.;
                    if (fvind > 0) {
                        gap0 = lam - P.values[fvind-1];
                    }
                    if (fvind < n-1) {
                        gap1 = P.values[fvind+1] - lam;
                    }
                    mpi_world_printf("    Fiedler vector:\n");
                    mpi_world_printf("      ind=%i\n", fvind);
                    mpi_world_printf("      lambda=%.18g\n", lam);
                    mpi_world_printf("      gap0=%8.3e\n", gap0);
                    mpi_world_printf("      gap1=%8.3e\n", gap1);
                    if (root) {
                        write_data_safely("fiedler", "f", 2,
                            opts.fiedler_vector_filename.c_str(),  
                            &fiedler[0], n, 1, false);
                    }
                }
            }

            if (opts.residuals) {
                mpi_world_printf("    Computing residuals.\n");
                std::vector<double> resids;
                tlist.start_event("residuals");
                P.residuals(resids);
                tlist.end_event();
                if (root) { tlist.report_event(6); }
                
                if (root) {
                    write_data_safely("residuals", "r", 2,
                        opts.residuals_filename.c_str(), 
                        &resids[0], n, 1, false);
                }
                
                if (root) {
                    find_large_residuals(resids, P.values, 2.2e-16*P.n);
                }
            }
            
            if (opts.iparscores) {
                std::vector<double> ipars;
                tlist.start_event("iparscores");
                P.Z.inverse_participation_ratios(ipars, true);  
                tlist.end_event();
                if (root) { tlist.report_event(6); }
                
                if (root) {
                    write_data_safely("iparscores", "p", 2,
                        opts.ipar_filename.c_str(), &ipars[0], n, 1, false);
                }
            }
            
            if (opts.matrix == opts.laplacian_matrix && opts.commute) {
                compute_commute_times(A, P);
            }
        }
        
        // handle Markov matrix
        if (opts.matrix == opts.normalized_laplacian_matrix && opts.markov) {
            output_markov_data(g, A, P);
        }
        
        // 
    }
    
    return 0;
} 

int main_blacs(int argc, char **argv, int nprow, int npcol) 
{
    int ictxt;
    int myrow, mycol;
    
    blacs_init();
    
    Cblacs_get(-1,0,&ictxt);
    Cblacs_gridinit(&ictxt, "R", npcol, npcol);    
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol );
    
    if (myrow == -1) {
        return (0);
    }
    
    bool root = (myrow == 0 && mycol == 0);
    
    // check that this process is in the computational grid and assert otherwise.
    assert(myrow < nprow);
    assert(mycol < npcol);
    
    triplet_data g;
    
    if (root) {
        const char* gfile = opts.graph_filename.c_str();
        printf("Loading graph %s.\n", gfile);
        tlist.start_event("load_graph");
        if (g.read_smat(gfile) == false) {
            printf("Error reading %s.\n", gfile);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        tlist.end_event(); tlist.report_event(2);
        
        printf("\n");
        printf("Graph information\n");
        printf("  Vertices: %i\n", g.nrows);
        printf("  Edges: %i\n", g.nnz);
        printf("\n");
        
        if (opts.verbose) {
            printf("Checking graph data.\n");
        }
        if (g.min_degree() < 1) {
            printf("Error: the graph has a disconnected vertex.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (!g.is_symmetric()) {
            printf("Error: the graph is not-symmetric.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (!g.has_repeated_edges()) {
            printf("Error: the graph has repeated edges.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        if (opts.verbose) {
            printf("Broadcasting graph data.\n");
        }
        
    }
    
    broadcast_triplet_data(ictxt, g);
    
    // we now allocate memory for the dense matrix to store the graph
    scalapack_distributed_matrix A;
    int n = g.nrows;
    if (opts.nb > n) {
        int newnb = opts.nb;
        
        while (newnb>1 && newnb>n/nprow) {
            newnb = newnb/2;
        }
        if (opts.verbose && root) {
            printf("Adjusting nb=%i to nb=%i because n=%i\n", 
                opts.nb, newnb, n);
        }
        opts.nb = newnb;
    }
    A.init(ictxt, n, n, opts.nb, opts.nb);
    
    tlist.start_event("construct_matrix");
    
    // allocate local storage for the matrix
    if (opts.verbose) {
        printf("[%3i x %3i] allocating %Zu bytes for A (n=%i, mb=%i, nb=%i, ap=%i, aq=%i)\n",
            myrow, mycol, A.bytes(),
            A.n, A.mb, A.nb, A.ap, A.aq );
    }
        
    if (!A.allocate()) {
        printf("[%3i x %3i] couldn't allocate memory, exiting...", 
           myrow, mycol);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    assign_graph_matrix(g, A, opts.matrix);
    
    tlist.end_event();
    if (root) {
        tlist.report_event(2);
    }
    if (root) { printf("\n"); }
        
    main_compute<scalapack_symmetric_eigen>(root, g, A);
    
    if (root) {
        printf("\n");
        printf("Done.\n");
        printf("  total_compute time: %.1f sec.\n", tlist.ttime);
        printf("\n");
    }
    
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
