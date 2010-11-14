
#ifdef __cplusplus
extern "C" {
#endif    


extern void   pdelset_( double *A, int *ia, int *ja, int *desca, double *alpha);
extern double pdlamch_( int *ictxt, char *cmach);
extern int    indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern int    indxg2l_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern int    indxl2g_( int *indxloc, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern int    numroc_( int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern void   descinit_( int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc,
                                int *ictxt, int *lld, int *info);
extern void   pdlaset_( char *uplo, int *m, int *n, double *alpha, double *beta, double *A, int *ia, int *ja, int *descA );
extern double pdlange_( char *norm, int *m, int *n, double *A, int *ia, int *ja, int *desca, double *work);
extern void   pdlacpy_( char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca,
                                double *b, int *ib, int *jb, int *descb);
extern void   pdgesv_( int *n, int *nrhs, double *A, int *ia, int *ja, int *desca, int* ipiv,
                                double *B, int *ib, int *jb, int *descb, int *info);
extern void   pdgesvd_( char *jobu, char *jobvt, int *m, int *n, double *a, int *ia, int *ja, int *desca,
                                double *s, double *u, int *iu, int *ju, int *descu,
                                double *vt, int *ivt, int *jvt, int *descvt, double *work, int *lwork, int *info);
extern void   pdgemm_( char *TRANSA, char *TRANSB, int * M, int * N, int * K, double * ALPHA,
                                double * A, int * IA, int * JA, int * DESCA, double * B, int * IB, int * JB, int * DESCB,
                                double * BETA, double * C, int * IC, int * JC, int * DESCC );
extern int    indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);

extern void   infog2l_( int *grind, int *gcind, int *desc, int *nprow, int *npcol,
                        int *myrow, int *mycol, int *lrind, int *lcind, int *rsrc, int *csrc);

extern void pdsymm_( char* side, char* uplo, int *M, int *N, double *ALPHA, 
    double *A, int *IA, int *JA, int *DESCA,
    double *B, int *IB, int *JB, int *DESCB, 
    double *BETA, 
    double *C, int *IC, int *JC, int *DESCC );


extern void pdsyev_( char *jobz, char *uplo, int *n,
                double *a, int *ia, int *ja, int *desca, double *w, double *z, int *iz, int *jz, int *descz,
                double *work, int *lwork, int *info );

extern void dcombnrm2_(double *x, double *y);
extern void pdtreecomb_(int *ictxt, char *scope, int *n, double *mine, int *rdest0, int *cdest0, void (*func)(double *x, double *y));

extern void pdlaprnt_(int *m, int *n, double *a, int *ia, int *ja, int *desca,
    int *iprnt, int *jprnt, char *cname,  int *nout, double *work, int clen);
    
extern void pdlawrite_( const char* filenam, int *m, int *n, double *a,
                        int *ia, int *ja, int *desca, 
                        int *piwrite, int *pjwrite, double *work,
                        int filenamelen );    
    
extern void pdlared2d_(int *m, int *ia, int *ja, int *desc, double *byrow, double *byall, double *work, int *lwork);
extern void pdlared1d_(int *n, int *ia, int *ja, int *desc, double *bycol, double *byall, double *work, int *lwork);


int sl_gridreshape_(int *ctxt, int *pstart, int *row_major_in, int *row_major_out, int *P, int *Q);

#ifdef __cplusplus
};
#endif
