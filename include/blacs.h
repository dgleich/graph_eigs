
#ifdef __cplusplus
extern "C" {
#endif    

extern void   Cblacs_pinfo( int* mypnum, int* nprocs);
extern void   Cblacs_get( int context, int request, int* value);
extern int    Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
extern void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
extern void   Cblacs_gridexit( int context);
extern void   Cblacs_exit( int error_code);
extern void   Cblacs_barrier(int ConTxt, char *scope);
extern void   Cigebs2d(int ConTxt, char *scope, char *top, int m, int n, int *A, int lda);
extern void   Cdgebs2d(int ConTxt, char *scope, char *top, int m, int n, double *A, int lda);
extern void   Cigebr2d(int ConTxt, char *scope, char *top, int m, int n, int *A,
              int lda, int rsrc, int csrc);
extern void   Cdgebr2d(int ConTxt, char *scope, char *top, int m, int n, double *A,
              int lda, int rsrc, int csrc);



extern void blacs_pinfo_( int* mypnum, int* nprocs);
extern void blacs_get_(int *context, int *what, int *val);
extern void blacs_gridinit_(int *context, char *order, int *nprow, int*npcol);
extern void blacs_gridexit_(int *context);
extern void blacs_gridinfo_( int *context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
extern void blacs_exit_( int* error_code);

#ifdef __cplusplus
};
#endif
