/**
 * @file scalapack_symmetric_eigen.cc
 * Compute eigenvalues and eigenvectors of symmetric problems with scalapack
 * @author David F. Gleich
 * 
 * Notes
 * -----
 * 
 * This code is rapidly evolving and the interface and abstractions haven't
 * settled down yet.  Thus, documentation will be scattershot until that
 * happens.
 * 
 * Todo
 * ----
 * TODO add debug mode for small problems to display intermediate data
 * TODO check that all calloced memory is freed
 * 
 * History
 * -------
 * :2010-11-08: Initial version
 * :2010-11-09: Added residuals
 * :2010-11-10: Added initial tridiagonal reduction
 */

#include <blas.h>
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
                     
extern void pdsyevr_tri_( char* jobc, char *jobz, char *range, char *uplo, int *n, 
                     double *a, int *ia, int *ja, int *desca, 
                     double *vl, double *vu, int *il, int *iu, int *m, int *nz,
                     double *w, double *z, int *iz, int *jz, int *descz, 
                     double *work, int *lwork, int *iwork, int *liwork,
                     int *jstate, int *info );  
                     
extern void pdsyev_tri_( char* jobc, char *jobz, char *uplo, int *n, 
                     double *a, int *ia, int *ja, int *desca, 
                     double *w, double *z, int *iz, int *jz, int *descz, 
                     double *work, int *lwork,
                     int *jstate, int *info );                       
                     
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
        ap = numroc_(&m, &mb, &myrow, &izero, &nprow);
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
    
    void print(int pi, int pj, char* cname, size_t cnamemax) {
        int izero = 0, ione = 1, isix = 6;
        std::vector<double> work(aq);
        
        int clen = strnlen(cname, cnamemax);
        
        pdlaprnt_(&m, &n, A, &ione, &ione, desc, 
            &izero, &izero, cname, &isix, &work[0], clen);
    }
    
    void write(std::string filename, int pi, int pj) {
        int izero = 0, ione = 1;
        std::vector<double> work(mb);
        const char *fnstr = filename.c_str();
        int fnlen = filename.size();
        pdlawrite_(fnstr, &m, &n, A, &ione, &ione, desc,
            &izero, &izero, &work[0], fnlen);
    }
    
    /** All processors set all of their entries to a constant. */
    void set_to_constant(double a) {
        off_t nel = (off_t)ap*(off_t)aq;
        for (off_t i=0; i<nel; ++i) {
            A[i] = a;
        }
    }

    void local2global(int li, int lj, int& gi, int& gj) {
        int izero=0;
        li+=1; // convert to fortran indices
        lj+=1;
        gi = indxl2g_(&li, &mb, &myrow, &izero, &nprow) - 1;
        gj = indxl2g_(&lj, &nb, &mycol, &izero, &npcol) - 1;
    }
    
    /**
     * @param gi the global row index
     * @param gj the global column index
     * @param [output] li the local row index
     * @param [output] lj the local column index
     * @return true if the element is local, or false otherwise.
     */
    bool global2local(int gi, int gj, int &li, int &lj) {
        int iarow, iacol;
        gi += 1;// adjust to fortran indices
        gj += 1; // adjust to fortran indices
        
        infog2l_(&gi, &gj, desc, &nprow, &npcol, &myrow, &mycol, 
            &li, &lj, &iarow, &iacol);
        li -= 1;
        lj -= 1;
        if (iarow == myrow && iacol == mycol) {
            return true;
        } else {
            return false;
        }        
    }
    
    /** Set A[gi,gj] = a.  
     * This command will only take effect on the processor that owns
     * gi,gj.
     */
    void set(int gi, int gj, double a) {
        gi += 1;// adjust to fortran indices
        gj += 1; // adjust to fortran indices
        pdelset_(A, &gi, &gj, desc, &a);
    }
    
    /** Set A[gi,gj] = A[gi,gj] + a
     * This command will only take effect on the processor 
     * that owns gi, gj.  It does not involve any communication.
     */
    void incr(int gi, int gj, double a) {
        int li, lj;
        if (global2local(gi, gj, li, lj)) {
            A[li+lj*ap] += a;
        }
    }
    
    
    /** Compute the norm of each column and aggregate on the first grid row.
     * @param norms 
     *   the output vector of column norms, only specified on the first
     *   row of the grid.
     * @param broadcast 
     *   if true, then broadcast the column norms to all processors
     *   if false, then norms is only modified on processor 1
     * 
     * This routine needs 2*aq work memory on each processor, and another 
     * n memory to hold the column norms on the first row.
     */
    void column_norms(std::vector<double>& norms) {
        int izero=0,ione=1;
        // TODO improve the algorithm here to use a more stable accumulation
        // TODO implement broadcast to send the results to all processors
        std::vector< double > colnorms(aq,0.);
        for (int j=0; j<aq; ++j) {
            //double scale=0., sumsq=0.;
            //dlassq_(&A.ap, &A.A[j*Aq], &ione, &scale, &sumsq);
            //colnorms[j] = scale*scale*sumsq;
            colnorms[j] = dnrm2_(&ap, &A[j*ap], &ione);
            //printf("[%3i x %3i] colnorms[%i] = %g\n", myrow, mycol, j, colnorms[j]);
            colnorms[j] *= colnorms[j];
        }
        
        
        Cdgsum2d(ictxt, "Columnwise", " ", 
          1, aq, &colnorms[0], 1, 0, mycol);
    
        
        // we can't use treecomb here because it doesn't work with vectors
        // pdtreecomb_(&ictxt, "Columnwise", 
        //    &aq, &colnorms[0], &izero, &mycol,
        //    &dcombnrm2_);
        
        //printf("postcomb\n");
        for (int j=0; j<aq; ++j) {
            colnorms[j] = sqrt(colnorms[j]);
            //printf("[%3i x %3i] colnorms[%i] = %g\n", myrow, mycol, j, colnorms[j]);
        }
        
        // TODO fix this code to only write to processor 0
        if (myrow == 0) {
            norms.resize(n);
            int lwork = numroc_(&n, &nb, &izero, &izero, &npcol);
            std::vector<double> work(lwork,0.);
            pdlared1d_( &n, &ione, &ione, desc, &colnorms[0], &norms[0],
                &work[0], &lwork);
            if (mycol == 0) {
                //printf("colnorms\n");
                for (int j=0; j<n; ++j) {
                    //printf("[%3i x %3i] norms[%i] = %g\n", myrow, mycol, j, norms[j]);
                }
            }
        }
    }
    
    /** Compute the inverse participation ratios for the column of a matrix.
     * 
     * For a vector, it's inverse participation score is:
     *  ipar(v) = sum(abs(v).^4)/((sum(abs(v)).^2)^2)
     * and 1/ipar(v) is the number of effective non-zeros in the vector.
     * 
     * Consider uniform v, then ipar(v) = 1/n.  For v with a single
     * nonzero, then ipar(v) = 1.
     * 
     * If the vectors are already normalized -- for instance when they are 
     * eigenvectors, which are usually taken to be normalized --
     * then we do not need to compute the denominator and save the work
     * and memory in this case.
     * 
     * @param ipars 
     *   the output vector of inverse participation scores, only specified on the first
     *   row of the grid.
     * 
     * @param normalized
     *   if true, then the matrix should already have normalized columns
     * 
     * This routine needs 2*aq work memory on each processor, and another 
     * n memory to hold the par scores on the first processor row.
     * If normalized is true, then it needs another ap work memory
     * on each processor.
     */
    void inverse_participation_ratios(std::vector<double>& ipars, 
        bool normalized) {
        int izero=0.,ione=1;
        std::vector<double> colnorms;
        
        if (normalized == false) {
            colnorms.resize(aq);
            for (int j=0; j<aq; ++j) {
                colnorms[j] = dnrm2_(&ap, &A[j*ap], &ione);
                colnorms[j] *= colnorms[j];
            }
            Cdgsum2d(ictxt, "Columnwise", " ", 
              1, aq, &colnorms[0], 1, 0, mycol);

        }

        // TODO implement broadcast to send the results to all processors
        std::vector< double > localpar(aq,0.);
        for (int j=0; j<aq; ++j) {
            for (int i=0; i<ap; ++i) {
                double score = fabs(A[j*ap+i]);
                score *= score;
                score *= score;
                localpar[j] += score;
            }
            //printf("[%3i x %3i] localpar[%i] = %g\n", myrow, mycol, j, localpar[j]);
        }

        Cdgsum2d(ictxt, "Columnwise", " ", 
          1, aq, &localpar[0], 1, 0, mycol);
        
        // TODO fix this code to only write to processor 0
        if (myrow == 0) {
            // rescale the normalized scores
            if (normalized == false) {
                for (int j=0; j<aq; ++j) {
                    localpar[j] /= (colnorms[j]*colnorms[j]);
                }
            }
            //printf("postcomb\n");
            for (int j=0; j<aq; ++j) {
                //printf("[%3i x %3i] localpar[%i] = %g\n", myrow, mycol, j, localpar[j]);
            }
            ipars.resize(n);
            int lwork = numroc_(&n, &nb, &izero, &izero, &npcol);
            std::vector<double> work(lwork,0.);
            pdlared1d_( &n, &ione, &ione, desc, &localpar[0], &ipars[0],
                &work[0], &lwork);
            if (mycol == 0) {
                //printf("ipars\n");
                for (int j=0; j<m; ++j) {
                    //printf("[%3i x %3i] ipars[%i] = %g\n", myrow, mycol, j, ipars[j]);
                }
            }
        }
    }
    
    /** Scale the column of a distributed matrix.
     * For each column A(:,i), this function sets
     *   A(:,i) <- alpha*A(:,i)*s(i)
     * where s(i) is the ith entry in the column scaling
     * vector.  This vector must be available on all
     * processors currently.
     * 
     * @param s the vector of scaling factors.  s must be an array
     * of length n (the number of columns) and must be available
     * on all processors.
     * @param alpha an additional scaling factor.
     */
    void scale_columns(double *s, double alpha) {
        int gi, gj;
        for (int i=0; i<ap; ++i) {
            for (int j=0; j<aq; ++j) {
                local2global(i, j, gi, gj);
                A[i+j*ap] *= alpha*s[gj];
            }
        }
    }
    
    /** Scale the rows of a distributed matrix.
     * For each row A(i,:), this function sets
     *   A(i,:) <- alpha*A(i,:)*s(i)
     * where s(i) is the ith entry in the row scaling
     * vector.  This vector must be available on all
     * processors currently.
     * 
     * @param s the vector of scaling factors.  s must be an array
     * of length m (the number of rows) and must be available
     * on all processors.
     * @param alpha an additional scaling factor.
     */
    void scale_rows(double *s, double alpha) {
        int gi, gj;
        for (int i=0; i<ap; ++i) {
            for (int j=0; j<aq; ++j) {
                local2global(i, j, gi, gj);
                A[i+j*ap] *= alpha*s[gi];
            }
        }
    }
};

struct scalapack_symmetric_eigen {
    int ictxt;
    int n;
    scalapack_distributed_matrix A;
    double *work;
    int lwork;
    int *iwork;
    int liwork;
    int jstate[3];
    double tau1;
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
        _set_worksize();
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
            double worksize[5];
            int iworksize[5];
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
                
            lwork = (int)worksize[0]+4;
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
            return 5 + 5*n + nextra + (2 + nz)*nn;
        } else {            
            int nextra = 12*nn;
            if ((nb * (np0 + 1)) > nextra) {
                nextra = (nb * (np0 + 1));
            }
            return 5 + 5*n + nextra;
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
    
    void tridiag_reduce() {
        // TODO add code to check for valid setup
        int ione=1, izero=0, nvals=0, nvecs=0, info=-1;
        double dzero=0.;
        char* jobz="N";
        if (vectors) {
            jobz="V";
        }
        
        // save diagonal
        _save_diagonal();
        
        pdsyevr_tri_("T",jobz,"A","U",
            &n, A.A, &ione, &ione, A.desc, // A
            &dzero, &dzero, &izero, &izero, // eigenvalue range, not specified
            &nvals, &nvecs, // eigenvalue output size
            values, // eigenvalue output
            Z.A, &ione, &ione, Z.desc, // eigenvector output
            work, &lwork, iwork, &liwork, // workspace
            jstate, // job state
            &info);
            
        assert(info == 0);
            
        if (tridiag) {
            // save the tridiagonal factors
            for (int i=0; i<n; ++i) {
                T[i] = work[jstate[0]+i-1];
                T[i+n] = work[jstate[1]+i+jstate[2]-1];
            }
            T[2*n-1] = 0.; // set the final entry to 0
        }
        
    }
    
    void tridiag_compute() {
        
        int ione=1, izero=0, nvals=0, nvecs=0, info=-1;
        double dzero=0.;
        char* jobz="N";
        if (vectors) {
            jobz="V";
        }
        
        pdsyevr_tri_("D",jobz,"A","U",
            &n, A.A, &ione, &ione, A.desc, // A
            &dzero, &dzero, &izero, &izero, // eigenvalue range, not specified
            &nvals, &nvecs, // eigenvalue output size
            values, // eigenvalue output
            Z.A, &ione, &ione, Z.desc, // eigenvector output
            work, &lwork, iwork, &liwork, // workspace
            jstate, // job state
            &info);
            
        _restore_diagonal();
            
        assert(info == 0);
        assert(nvals == n);
        if (vectors) {
            assert(nvecs == n);
        }
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
 
    /** Compute the residuals for the eigenvectors 
     * 
     * This computes
     *   residnorms(i) = ||A*v - v*l||
     * 
     * This will only store the vector of residuals on the 
     * first row of the processor grid.
     * 
     * This function does not require the matrix to be 
     * reloaded if both halves (upper and lower) were initially
     * specified.  
     */
    void residuals(std::vector<double>& residnorms) {
        double done = 1.;
        int ione=1;
        // repurpose lwork
        assert(lwork > A.ap*A.aq);
        
        scalapack_distributed_matrix resids;
        resids.init(ictxt, n, n, A.nb, A.nb);
        resids.A = work;
        
        // compute -V*L
        memcpy(work, Z.A, sizeof(double)*Z.ap*Z.aq);
        resids.scale_columns(values, -1.0);
        
        // compute A*V-V*L for all the eigenvectors
        pdsymm_("L", "L", 
            &n, &n, &done,
            A.A, &ione, &ione, A.desc,
            Z.A, &ione, &ione, Z.desc,
            &done,
            work, &ione, &ione, Z.desc);
            
        
        resids.column_norms(residnorms);
    }
    
    /** Compute the residuals for the eigenvectors 
     * 
     * This computes
     *   residnorms(i) = ||A*v - v*l||
     * 
     * This will only store the vector of residuals on the 
     * first row of the processor grid.
     * 
     * This function requires the matrix A to be reloaded.  
     * It uses pdgemm instead of pdsymm.
     */
    void reloaded_residuals(std::vector<double>& residnorms) {
        double done = 1.;
        int ione=1;
        // repurpose lwork
        assert(lwork > A.ap*A.aq);
        
        scalapack_distributed_matrix resids;
        resids.init(ictxt, n, n, A.nb, A.nb);
        resids.A = work;
        
        // compute -V*L
        memcpy(work, Z.A, sizeof(double)*Z.ap*Z.aq);
        resids.scale_columns(values, -1.0);
        
        // compute A*V-V*L for all the eigenvectors
        pdgemm_("N", "N", 
            &n, &n, &n, &done,
            A.A, &ione, &ione, A.desc,
            Z.A, &ione, &ione, Z.desc,
            &done,
            work, &ione, &ione, Z.desc);
            
        
        resids.column_norms(residnorms);
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
 
class scalapack_symmetric_eigen_std {
public:    
    int ictxt;
    int n;
    scalapack_distributed_matrix A;
    double *work;
    int lwork;
    int jstate[3];
    bool vectors;
    bool minmemory;
    bool tridiag;
    
    bool myZ;
    scalapack_distributed_matrix Z;
    
    double *T;
    double *values;
    double *diag;
    
    scalapack_symmetric_eigen_std(scalapack_distributed_matrix& A_)
        : ictxt(A_.ictxt), n(A_.n), A(A_),
          work(NULL), lwork(0),
          vectors(false), minmemory(false), tridiag(false),
          myZ(true),
          T(NULL), values(NULL), diag(NULL)
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
        _set_worksize();
    }
    
    /** Designed for sub-classes to override with any extra memory they need. */ 
    virtual size_t _bytes() {
        return 0;
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
        size_t nwork = sizeof(double)*(lwork + A.ap) + _bytes();
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
        
        return true;
    }
    
    virtual void _set_worksize() {
        // do a workspace query
        int ione=1, lworkq=-1, info;
        double worksize[5];
        char* jobz="N";
        if (vectors) {
            jobz="V";
        }
           
        pdsyev_(jobz, "U", 
            &n, 
            NULL, &ione, &ione, A.desc, // A
            NULL, // eigenvalue output
            NULL, &ione, &ione, A.desc, // Z (but not allocated)
            worksize, &lworkq, &info);
            
        lwork = (int)worksize[0]+3;
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
    
    int _nprocs() {
        int nprow, npcol, myrow, mycol;
        
        Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
        return nprow*npcol;
    }
    
    virtual void tridiag_reduce() {
        // TODO add code to check for valid setup
        int ione=1, info=-1;
        char* jobz="N";
        if (vectors) {
            jobz="V";
        }
        
        // save diagonal
        _save_diagonal();
        
        pdsyev_tri_("T",jobz,"U",
            &n, A.A, &ione, &ione, A.desc, // A
            values, // eigenvalue output
            Z.A, &ione, &ione, Z.desc, // eigenvector output
            work, &lwork,
            jstate, // job state
            &info);
            
        assert(info == 0);
            
        if (tridiag) {
            // save the tridiagonal factors
            for (int i=0; i<n; ++i) {
                T[i] = work[jstate[0]+i-1];
                T[i+n] = work[jstate[1]+i+jstate[2]-1];
            }
            T[2*n-1] = 0.; // set the final entry to 0
        }
        
    }
    
    virtual void tridiag_compute() {
        
        int ione=1,info=-1;
        char* jobz="N";
        if (vectors) {
            jobz="V";
        }
        
        pdsyev_tri_("D",jobz,"U",
            &n, A.A, &ione, &ione, A.desc, // A
            values, // eigenvalue output
            Z.A, &ione, &ione, Z.desc, // eigenvector output
            work, &lwork, 
            jstate, // job state
            &info);
            
        _restore_diagonal();
            
        assert(info == 0);
    }
    
    virtual void compute() {
        // TODO add code to check for valid setup
        
        int ione=1, info;
        char* jobz="N";
        if (vectors) {
            jobz="V";
        }
        
        
        
        // save diagonal
        _save_diagonal();
        
        pdsyev_(jobz,"U",
            &n, A.A, &ione, &ione, A.desc, // A
            values, // eigenvalue output
            Z.A, &ione, &ione, Z.desc, // eigenvector output
            work, &lwork, // workspace
            &info);
            
                        
        _restore_diagonal();
            
        assert(info == 0);
    }
 
    /** Compute the residuals for the eigenvectors 
     * 
     * This computes
     *   residnorms(i) = ||A*v - v*l||
     * 
     * This will only store the vector of residuals on the 
     * first row of the processor grid.
     * 
     * This function does not require the matrix to be 
     * reloaded if both halves (upper and lower) were initially
     * specified.  
     */
    void residuals(std::vector<double>& residnorms) {
        double done = 1.;
        int ione=1;
        // repurpose lwork
        assert(lwork > A.ap*A.aq);
        
        scalapack_distributed_matrix resids;
        resids.init(ictxt, n, n, A.nb, A.nb);
        resids.A = work;
        
        // compute -V*L
        memcpy(work, Z.A, sizeof(double)*Z.ap*Z.aq);
        resids.scale_columns(values, -1.0);
        
        // compute A*V-V*L for all the eigenvectors
        pdsymm_("L", "L", 
            &n, &n, &done,
            A.A, &ione, &ione, A.desc,
            Z.A, &ione, &ione, Z.desc,
            &done,
            work, &ione, &ione, Z.desc);
            
        
        resids.column_norms(residnorms);
    }
    
    /** Compute the residuals for the eigenvectors 
     * 
     * This computes
     *   residnorms(i) = ||A*v - v*l||
     * 
     * This will only store the vector of residuals on the 
     * first row of the processor grid.
     * 
     * This function requires the matrix A to be reloaded.  
     * It uses pdgemm instead of pdsymm.
     */
    void reloaded_residuals(std::vector<double>& residnorms) {
        double done = 1.;
        int ione=1;
        // repurpose lwork
        assert(lwork > A.ap*A.aq);
        
        scalapack_distributed_matrix resids;
        resids.init(ictxt, n, n, A.nb, A.nb);
        resids.A = work;
        
        // compute -V*L
        memcpy(work, Z.A, sizeof(double)*Z.ap*Z.aq);
        resids.scale_columns(values, -1.0);
        
        // compute A*V-V*L for all the eigenvectors
        pdgemm_("N", "N", 
            &n, &n, &n, &done,
            A.A, &ione, &ione, A.desc,
            Z.A, &ione, &ione, Z.desc,
            &done,
            work, &ione, &ione, Z.desc);
            
        
        resids.column_norms(residnorms);
    }
    
    /** Free any memory allocated in this operation. */
    virtual void free() {
        if (vectors && myZ) {
            Z.free();
        }
        if (tridiag) {
            ::free(T);
        }
        ::free(values);
        ::free(work);
    }
}; 

class scalapack_symmetric_eigen_mrrr : scalapack_symmetric_eigen_std {
public:    
    int *iwork;
    int liwork;
    
    scalapack_symmetric_eigen_mrrr(scalapack_distributed_matrix& A_)
        : scalapack_symmetric_eigen_std(A_), iwork(NULL), liwork(0)
    {
    }
    
    size_t _bytes() {
        return sizeof(int)*liwork;
    }
    
    bool allocate() {
        if (scalapack_symmetric_eigen_std::allocate()) {
            iwork = (int*)calloc(liwork, sizeof(int));
            return true;
        } else {
            return false;
        }
    }
    
    /** Free any memory allocated in this operation. */
    void free() {
        scalapack_symmetric_eigen_std::free();
        ::free(iwork);
    }
    
    void _set_worksize() {
        
        if (minmemory) {
            lwork = (int)_workmin();
            liwork = (int)_iworkmin();
        } else {
            // do a workspace query
            int ione=1, izero=0, nvals, nvecs, lworkq=-1, info;
            double worksize[5];
            int iworksize[5];
            double dzero=0.;
            char* jobz="N";
            if (vectors) {
                jobz="V";
            }
               
            pdsyevr_tri_("A", jobz, "A", "U", 
                &n, 
                NULL, &ione, &ione, A.desc, // A
                &dzero, &dzero, &izero, &izero, // eigenvalue range, not specified
                &nvals, &nvecs, // eigenvalue output size
                NULL, // eigenvalue output
                NULL, &ione, &ione, A.desc, // Z (but not allocated)
                worksize, &lworkq, iworksize, &lworkq,
                jstate,
                &info);
                
            lwork = (int)worksize[0];
            liwork = (int)iworksize[0];
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
            return 5 + 5*n + nextra + (2 + nz)*nn;
        } else {            
            int nextra = 12*nn;
            if ((nb * (np0 + 1)) > nextra) {
                nextra = (nb * (np0 + 1));
            }
            return 5 + 5*n + nextra;
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
    
    void tridiag_reduce() {
        // TODO add code to check for valid setup
        int ione=1, izero=0, nvals=0, nvecs=0, info=-1;
        double dzero=0.;
        char* jobz="N";
        if (vectors) {
            jobz="V";
        }
        
        // save diagonal
        _save_diagonal();
        
        pdsyevr_tri_("T",jobz,"A","U",
            &n, A.A, &ione, &ione, A.desc, // A
            &dzero, &dzero, &izero, &izero, // eigenvalue range, not specified
            &nvals, &nvecs, // eigenvalue output size
            values, // eigenvalue output
            Z.A, &ione, &ione, Z.desc, // eigenvector output
            work, &lwork, iwork, &liwork, // workspace
            jstate, // job state
            &info);
            
        assert(info == 0);
            
        if (tridiag) {
            // save the tridiagonal factors
            for (int i=0; i<n; ++i) {
                T[i] = work[jstate[0]+i-1];
                T[i+n] = work[jstate[1]+i+jstate[2]-1];
            }
            T[2*n-1] = 0.; // set the final entry to 0
        }
        
    }
    
    void tridiag_compute() {
        
        int ione=1, izero=0, nvals=0, nvecs=0, info=-1;
        double dzero=0.;
        char* jobz="N";
        if (vectors) {
            jobz="V";
        }
        
        pdsyevr_tri_("D",jobz,"A","U",
            &n, A.A, &ione, &ione, A.desc, // A
            &dzero, &dzero, &izero, &izero, // eigenvalue range, not specified
            &nvals, &nvecs, // eigenvalue output size
            values, // eigenvalue output
            Z.A, &ione, &ione, Z.desc, // eigenvector output
            work, &lwork, iwork, &liwork, // workspace
            jstate, // job state
            &info);
            
        _restore_diagonal();
            
        assert(info == 0);
        assert(nvals == n);
        if (vectors) {
            assert(nvecs == n);
        }
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
}; 
 
