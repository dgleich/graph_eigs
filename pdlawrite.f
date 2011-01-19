      SUBROUTINE PDLAWRITE( FILENAM, M, N, A, IA, JA, DESCA, IRWRIT,
     $                      ICWRIT, WORK )
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     August 12, 2001
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    FILENAM
      INTEGER            IA, ICWRIT, IRWRIT, JA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * ), WORK( * )
*     ..
*
* Purpose
* =======
*
* PDLAWRITE writes to a file named FILNAMa distributed matrix sub( A )
* denoting A(IA:IA+M-1,JA:JA+N-1). The local pieces are sent to and
* written by the process of coordinates (IRWWRITE, ICWRIT).
*
* WORK must be of size >= MB_ = DESCA( MB_ ).
*
* Further Details
* ===============
*
* Contributed by Song Jin, University of Tennessee, 1996.
*
* =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            NOUT
      PARAMETER          ( NOUT = 13 )
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ISIOPROCESSOR
      INTEGER            CSRC, I, ICTXT, IEND, ISIZE, ISTART, J, JEND,
     $                   JSIZE, JSTART, LDD, LWORK, MB, MM, MYCOL,
     $                   MYROW, NB, NN, NPCOL, NPROW, RSRC
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Local Arrays ..
      INTEGER            DESCWORK( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, PDGEADD, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN
*     ..
*     .. Executable Statements ..
      LWORK = DESCA( MB_ )
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      ISIOPROCESSOR = ( ( MYROW.EQ.IRWRIT ) .AND. ( MYCOL.EQ.ICWRIT ) )
*
      MM = MAX( 1, MIN( M, LWORK ) )
      NN = MAX( 1, INT( LWORK / MM ) )
      MB = MM
      NB = NN
      RSRC = IRWRIT
      CSRC = ICWRIT
      LDD = MAX( 1, MM )
      CALL DESCSET( DESCWORK, MM, NN, MB, NB, RSRC, CSRC, ICTXT, LDD )
      IF( ISIOPROCESSOR ) THEN
         OPEN( NOUT, FILE = FILENAM, STATUS = 'UNKNOWN',
     $       FORM = 'FORMATTED', ACCESS = 'SEQUENTIAL', ERR = 50 )
         REWIND ( NOUT )
         WRITE( NOUT, FMT = *, ERR = 50 )M, N
      END IF
      DO 40 JSTART = JA, JA + N - 1, NN
         JEND = MIN( JA+N-1, JSTART+NN-1 )
         JSIZE = JEND - JSTART + 1
         DO 30 ISTART = IA, IA + M - 1, MM
            IEND = MIN( IA+M-1, ISTART+MM-1 )
            ISIZE = IEND - ISTART + 1
            ALPHA = ONE
            BETA = ZERO
            CALL PDGEADD( 'NoTrans', ISIZE, JSIZE, ALPHA, A, ISTART,
     $                    JSTART, DESCA, BETA, WORK, 1, 1, DESCWORK )
            IF( ISIOPROCESSOR ) THEN
               DO 20 J = 1, JSIZE
                  DO 10 I = 1, ISIZE
                     WRITE( NOUT, FMT = *, ERR = 50 )WORK( I+( J-1 )*
     $                  LDD )
   10             CONTINUE
   20          CONTINUE
            END IF
   30    CONTINUE
   40 CONTINUE
      IF( ISIOPROCESSOR ) THEN
         CLOSE ( NOUT, ERR = 50 )
      END IF
      WORK( 1 ) = DESCA( MB_ )
      RETURN
   50 CONTINUE
      CALL PXERBLA( DESCA( CTXT_ ), 'PLAWRITE', 1 )
      WORK( 1 ) = DESCA( MB_ )
      RETURN
*
*     End of PDLAWRITE
*
      END
