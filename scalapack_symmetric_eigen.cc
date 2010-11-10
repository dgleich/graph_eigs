/**
 * @file scalapack_symmetric_eigen.cc
 * Compute eigenvalues and eigenvectors of symmetric problems with scalapack
 * @author David F. Gleich
 */

#include <blacs.h>
#include <scalapack.h>

#ifdef __cplusplus
extern "C" {
#endif
extern void pdsyevr_( char *jobz, char *range, char *uplo, int *n, 
                     double *a, int *ia, int *ja, int *desca, 
                     double *vl, double *vu, int *il, int *iu, int *m, int *nz,
                     double *w, double *z, int *iz, int *jz, int *descz, 
                     double *work, int *lwork, int *iwork, int *liwork,
                     int *info );
#ifdef __cplusplus
};
#endif

struct scalapack_distributed_matrix {
    int ictxt;
    int desc[9];
    int m, n, mb, nb, ap, aq; // size and blocking
    int myrow, mycol, nprow, npcol;
    double *A;
    
    bool init(int ictxt_, int m_, int n_, int mb_, int nb_) {
        int izero = 0, info;
        ictxt = ictxt_;
        m = m_;
        n = n_;
        mb = mb_;
        nb = nb_;
        A = NULL;
        Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
        ap = numroc_(&n, &mb, &myrow, &izero, &nprow);
        aq = numroc_(&n, &nb, &mycol, &izero, &npcol);
        descinit_(desc, &m, &n, &mb, &nb, &izero, &izero, 
            &ictxt, &ap, &info);
        if (info != 0) {
            return false;
        } else {
            return true;
        }
    }

    
    size_t bytes() {
        return sizeof(double)*ap*aq;
    }
    
    bool allocate() {
        assert(A == NULL);
        A = (double*)calloc(ap*aq, sizeof(double));
        if (A) {
            return true;
        } else {
            return false;
        }
    }
    
    void free() {
        if (A) {
            ::free(A);
            A = NULL;
        }
    }

    void local2global(int li, int lj, int& gi, int& gj) {
        int izero=0;
        li+=1; // convert to fortran indices
        lj+=1;
        gi = indxl2g_(&li, &mb, &myrow, &izero, &nprow) - 1;
        gj = indxl2g_(&lj, &nb, &mycol, &izero, &npcol) - 1;
    }
};
 
/** Usage:
 * scalapack_distributed_matrix A(n, n, nb, nb);
 * A.allocate();
 * // assign(A);
 * scalapack_symmetric_eigen P(A);
 * P.setup(false // vectors, true // minimial memory, true // tridiagonal);
 * size_t memreq = P.bytes();
 * P.allocate();
 * P.tridiag(); // only does something if tridiag is true
 * P.compute();
 * P.get_values();
 * P.get_vectors();
 * P.get_tridiag();
 * 
 * Notes:
 * This class is designed to share as much information as possible
 *   so it won't copy the matrix A, and will modify it in place.
 * 
 */
struct scalapack_symmetric_eigen {
    int ictxt;
    int n;
    scalapack_distributed_matrix A;
    double *work;
    int lwork;
    int *iwork;
    int liwork;
    bool vectors;
    bool minmemory;
    bool tridiag;
    
    bool myZ;
    scalapack_distributed_matrix Z;
    
    double *T;
    double *values;
    double *diag;
    
    scalapack_symmetric_eigen(scalapack_distributed_matrix& A_)
        : ictxt(A_.ictxt), n(A_.n), A(A_),
          work(NULL), lwork(0), iwork(NULL), liwork(0),
          vectors(false), minmemory(false), tridiag(false),
          myZ(true),
          T(NULL), values(NULL)
    {
    }
    
    /** Call this function to manage the eigenvector matrix yourself. 
     * This function must be called before setup.
     */
    void set_eigenvector_matrix(scalapack_distributed_matrix Z_) {
        Z = Z_;
        myZ = false;
    }
    
    /** Pick which type of problem to solve.
     * If you want to use your own matrix for eigenvectors,
     * it must have already been set at this point.
     */
    void setup(bool vectors_, bool minmemory_, bool tridiag_) {
        vectors = vectors_;
        minmemory = minmemory_;
        tridiag = tridiag_;
        
        if (vectors && myZ) {
            Z.init(ictxt, n, n, A.nb, A.nb);
        }
    }
    
    size_t bytes() {
        size_t nvalues = n*sizeof(double);
        size_t nvectors = 0;
        if (vectors) {
            if (myZ) {
                nvectors = Z.bytes();
            }
        }
        _set_worksize();
        // we need to allocate an extra A.ap elements to store the diagonal
        // of A.
        size_t nwork = sizeof(double)*(lwork + A.ap) + sizeof(int)*liwork;
        size_t ntridiag = sizeof(double)*2*n;
        
        return nvalues + nvectors + nwork + ntridiag;
    }
    
    bool allocate() {
        if (vectors && myZ) {
            Z.allocate();
        }
        if (tridiag) {
            T = (double*)calloc(2*n, sizeof(double));
        }
        values = (double*)calloc(n, sizeof(double));
        diag = (double*)calloc(A.ap, sizeof(double));
        
        _set_worksize();
        work = (double*)calloc(lwork, sizeof(double));
        iwork = (int*)calloc(liwork, sizeof(int));
        
        return true;
    }
    
    void _set_worksize() {
        
        if (minmemory) {
            lwork = (int)_workmin();
            liwork = (int)_iworkmin();
        } else {
            // do a workspace query
            int ione=1, izero=0, nvals, nvecs, lworkq=-1, info;
            double worksize[3];
            int iworksize[3];
            double dzero=0.;
            char* jobz="N";
            if (vectors) {
                jobz="V";
            }
               
            pdsyevr_(jobz, "A", "U", 
                &n, 
                NULL, &ione, &ione, A.desc, // A
                &dzero, &dzero, &izero, &izero, // eigenvalue range, not specified
                &nvals, &nvecs, // eigenvalue output size
                NULL, // eigenvalue output
                NULL, &ione, &ione, A.desc, // Z (but not allocated)
                worksize, &lworkq, iworksize, &lworkq,
                &info);
                
            lwork = (int)worksize[0];
            liwork = (int)iworksize[0];
        }
    }
    
    void _save_diagonal() {
        int gi, gj;
        for (int i=0; i<A.ap; ++i) {
            for (int j=0; j<A.aq; ++j) {
                A.local2global(i, j, gi, gj);
                if (gi==gj) {
                    diag[i] = A.A[i+j*A.ap];
                }
            }
        }
    }
    
    void _restore_diagonal() {
        int gi, gj;
        for (int i=0; i<A.ap; ++i) {
            for (int j=0; j<A.aq; ++j) {
                A.local2global(i, j, gi, gj);
                if (gi==gj) {
                    A.A[i+j*A.ap] = diag[i];
                }
            }
        }
    }
    
    size_t _workmin() {
    
        int nprow, npcol, myrow, mycol;
        int izero=0;
        
        Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
      
        int nb = A.nb;
        
        int nn = 2;
        if (nb > nn) { nn = nb; }
        if (n > nn) { nn = n; }
        
        int np0 = numroc_( &nn, &nb, &izero, &izero, &nprow );
        int mq0 = numroc_( &nn, &nb, &izero, &izero, &npcol );
        
        if (vectors) {
            int nextra = np0*mq0 + 2*nb*nb;
            if (18*nn > nextra) {
                nextra = 18*nn; 
            }
            int nz = n / (nprow*npcol);
            if (n % (nprow*npcol) != 0) {
                nz += 1;
            }
            return 2 + 5*n + nextra + (2 + nz)*nn;
        } else {            
            int nextra = 12*nn;
            if ((nb * (np0 + 1)) > nextra) {
                nextra = (nb * (np0 + 1));
            }
            return 2 + 5*n + nextra;
        }
    }
    
    size_t _iworkmin() {
        int nnp = _nprocs()+1;
        if (nnp < 4) { nnp = 4; }
        if (nnp < n) { nnp = n; }
        
        if (vectors) {
            return (size_t)(12*nnp + 2*n);
        } else {
            return (size_t)(10*nnp + 2*n);
        }
    }
    
    int _nprocs() {
        int nprow, npcol, myrow, mycol;
        
        Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
        return nprow*npcol;
    }
    
    void compute() {
        // TODO add code to check for valid setup
        
        int ione=1, izero=0, nvals, nvecs, info;
        double dzero=0.;
        char* jobz="N";
        if (vectors) {
            jobz="V";
        }
        
        // save diagonal
        _save_diagonal();
        
        pdsyevr_(jobz,"A","U",
            &n, A.A, &ione, &ione, A.desc, // A
            &dzero, &dzero, &izero, &izero, // eigenvalue range, not specified
            &nvals, &nvecs, // eigenvalue output size
            values, // eigenvalue output
            Z.A, &ione, &ione, Z.desc, // eigenvector output
            work, &lwork, iwork, &liwork, // workspace
            &info);
            
                        
        _restore_diagonal();
            
        assert(info == 0);
        assert(nvals == n);
        if (vectors) {
            assert(nvecs == n);
        }
    }
    
    void _scale_eigenvectors_in_work() {
        int gi, gj;
        for (int i=0; i<Z.ap; ++i) {
            for (int j=0; j<Z.aq; ++j) {
                Z.local2global(i, j, gi, gj);
                    work[i+j*A.ap] *= -1.0*values[gj];
            }
        }
    }

    /** Compute the residuals for the eigenvectors */
    void residuals() {
        double done = 1., dzero = 0.;
        int ione=1;
        // repurpose lwork
        assert(lwork > A.ap*A.aq);
        //memset(work, 0, sizeof(double)*A.ap*A.aq);
        
        // compute -V*L
        memcpy(work, Z.A, sizeof(double)*Z.ap*A.aq);
        _scale_eigenvectors_in_work();
        
        // compute A*V-V*L for all the eigenvectors
        pdsymm_("L", "L", 
            &n, &n, &done,
            A.A, &ione, &ione, A.desc,
            Z.A, &ione, &ione, Z.desc,
            &done,
            work, &ione, &ione, Z.desc);
            
        // now compute norms of each column
        double maxresid = pdlange_("M", &n, &n, work, &ione, &ione, A.desc, &dzero);
        printf("Residual: %g\n", maxresid);
        
    }
    
    /** Free any memory allocated in this operation. */
    void free() {
        if (vectors && myZ) {
            Z.free();
        }
        if (tridiag) {
            ::free(T);
        }
        ::free(values);
        ::free(work);
        ::free(iwork);
    }
}; 
 
