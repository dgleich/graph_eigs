/**
 * @file scalapack_symmetric_eigen.cc
 * Compute eigenvalues and eigenvectors of symmetric problems with scalapack
 * @author David F. Gleich
 */

struct scalapack_distributed_matrix {
    int ictxt;
    int desc[9];
    int m, n, mb, nb, ap, aq; // size and blocking
    double *A;
    
    bool init(int ictxt_, int m_, int n_, int mb_, int nb_) {
        int izero = 0, info;
        ictxt = ictxt_;
        m = m_;
        n = n_;
        mb = mb_;
        nb = nb_;
        ap = numroc_(&n, &nblock, &myrow, &izero, &nprow);
        aq = numroc_(&n, &nblock, &mycol, &izero, &npcol);
        descinit_(desc, &m, &n, &mblock, &nblock, &izero, &izero, 
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
            free(A);
            A = NULL;
        }
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
    
    void scalapack_symmetric_eigen(scalapack_distributed_matrix& A_)
        : A(A_), ictxt(A_.ictxt), n(A_.n), 
          work(NULL), lwork(0), iwork(NULL), liwork(0),
          vectors(false), minmemory(false), tridiag(false),
          myZ(true),
          T(NULL), values(NULL);
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
            Z.init(ictxt, n, n, A_.nb, A_.nb);
        }
    }
    
    size_t bytes() {
        size_t nvalues = n*sizeof(double);
        size_t nvectors = 0;
        if (vectors) {
            if (myZ) {
                size_t nvectors = Z.bytes();
            }
        }
        _set_worksize();
        size_t nwork = sizeof(double)*lwork + sizeof(int)*liwork;
        size_t ntridiag = sizeof(double)*2*n;
        
        return nvalues + nvectors + nwork + ntridiag;
    }
    
    bool allocate() {
        if (vectors && myZ) {
            Z.allocate();
        }
        if (tridiag) {
            T = calloc(2*n, sizeof(double));
        }
        values = calloc(n, sizeof(double));
        
        _set_worksize();
        work = calloc(lwork, sizeof(double));
        iwork = calloc(liwork, sizeof(int));
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
    
    size_t _workmin() {
    
        int nprow, npcol, myrow, mycol;
        
        Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
      
        int nb = A.nb;
        
        int nn = 2;
        if (nb > nn) { nn = nb; }
        if (n > nn) { n = nn; }
        
        int np0 = numroc_( &nn, &nb, &izero, &izero, &nprow );
        int mp0 = numroc_( &nn, &nb, &izero, &izero, &npcol );
        
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
        int nnp = nprocs+1;
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

    /** Compute the residuals for the eigenvectors */
    void residuals() {}
    
    /** Free any memory allocated in this operation.
    void free() {
    }
}; 
 
