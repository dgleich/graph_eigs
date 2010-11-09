      SUBROUTINE PDSTEGR( JOBZ, RANGE, UPLO, N, D, E, 
     $                    VL, VU, IL, IU, M, NZ, W, Z, IZ,
     $                    JZ, DESCZ, WORK, LWORK, IWORK, LIWORK,
     $                    INFO )

      IMPLICIT NONE
*
*  -- ScaLAPACK routine (@(MODE)version *TBA*) --
*     University of California, Berkeley and
*     University of Tennessee, Knoxville. 
*     December 6, 2008
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, IZ, JZ, LIWORK, LWORK, M,
     $                   N, NZ
      DOUBLE PRECISION VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            DESCZ( * ), IWORK( * )
      DOUBLE PRECISION   D(*), E(*), W(*), WORK(*), Z(*)
*     ..
*
*  Purpose
*  =======
*
*  PDSTEGR computes selected eigenvalues and eigenvectors
*  of a real symmetric TRIDIAGONAL matrix by calling the recommended sequence
*  of ScaLAPACK routines.  Eigenvalues/vectors can be selected by
*  specifying a range of values or a range of indices for the desired
*  eigenvalues. 
*
*  THIS IS AN AUXILIARY CODE FOR STUDYING THE PARALLEL TRIDIAGONAL 
*  EIGENSOLVER!!!
*
*  Arguments
*  =========
*
*  JOBZ    (global input) CHARACTER*1
*          Specifies whether or not to compute the eigenvectors:
*          = 'V':  Compute eigenvalues and eigenvectors. 
*
*  RANGE   (global input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the interval [VL,VU] will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*  UPLO    (global input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns of the matrix A.  N >= 0
*
*  D,E     (local input) DOUBLE PRECISION arrays
*          both of local dimension (N), for the tridiagonal matrix.
*
*  VL      (global input) DOUBLE PRECISION 
*          If RANGE='V', the lower bound of the interval to be searched
*          for eigenvalues.  Not referenced if RANGE = 'A' or 'I'.
*
*  VU      (global input) DOUBLE PRECISION 
*          If RANGE='V', the upper bound of the interval to be searched
*          for eigenvalues.  Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (global input) INTEGER
*          If RANGE='I', the index (from smallest to largest) of the
*          smallest eigenvalue to be returned.  IL >= 1.
*          Not referenced if RANGE = 'A'.
*
*  IU      (global input) INTEGER
*          If RANGE='I', the index (from smallest to largest) of the
*          largest eigenvalue to be returned.  min(IL,N) <= IU <= N.
*          Not referenced if RANGE = 'A'.
*
*  M       (global output) INTEGER
*          Total number of eigenvalues found.  0 <= M <= N.
*
*  NZ      (global output) INTEGER
*          Total number of eigenvectors computed.  0 <= NZ <= M.
*          The number of columns of Z that are filled.
*          If JOBZ .EQ. 'V', NZ = M 
*
*  W       (global output) DOUBLE PRECISION array, dimension (N)
*          On normal exit, the first M entries contain the selected
*          eigenvalues in ascending order.
*
*  Z       (local output) DOUBLE PRECISION array,
*          global dimension (N, N),
*          local dimension ( LLD_Z, LOCc(JZ+N-1) )
*          If JOBZ = 'V', then on normal exit the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix
*          corresponding to the selected eigenvalues.
*
*  IZ      (global input) INTEGER
*          Z's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JZ      (global input) INTEGER
*          Z's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*          DESCZ( CTXT_ ) must equal DESCZ( CTXT_ )
*
*  WORK    (local workspace/output) DOUBLE PRECISION  array,
*          dimension (LWORK)
*          On return, WORK(1) contains the optimal amount of
*          workspace required for efficient execution.
*          If JOBZ='V' WORK(1) = optimal amount of workspace
*             required to compute eigenvalues and eigenvectors.
*
*  LWORK   (local input) INTEGER
*          Size of WORK
*          See below for definitions of variables used to define LWORK.
*          If eigenvectors are requested (JOBZ = 'V' ) then
*             the amount of workspace required is:
*             LWORK >= 2*N + 18*N +
*               (2 + ICEIL( NEIG, NPROW*NPCOL))*N
*          Note that this work space is smaller than the one required 
*          for the dense driver PDSYEVR.
*
*          Variable definitions:
*             NEIG = number of eigenvectors requested
*             ICEIL( X, Y ) is a ScaLAPACK function returning
*             ceiling(X/Y)
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the size
*          required for optimal performance for all work arrays. Each of
*          these values is returned in the first entry of the
*          corresponding work arrays, and no error message is issued by
*          PXERBLA.
*
*  IWORK   (local workspace) INTEGER array
*          On return, IWORK(1) contains the amount of integer workspace
*          required.
*
*  LIWORK  (local input) INTEGER
*          size of IWORK
*
*          Let  NNP = MAX( N, NPROW*NPCOL + 1, 4 ). Then:
*          LIWORK >= 12*NNP + 2*N when the eigenvectors are desired
*          LIWORK >= 10*NNP + 2*N when only the eigenvalues have to be computed
*          
*          If LIWORK = -1, then LIWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*     .. Parameters ..
      INTEGER            CTXT_
      PARAMETER          ( CTXT_ = 2 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, COLBRT, DOBCST, FINISH, FIRST, INDEIG,
     $                   LOWER, LQUERY, VALEIG, VSTART, WANTZ
      INTEGER            DOL, DOU, DSTCOL, DSTROW, EIGCNT, FRSTCL, I,
     $                   ICTXT, IIL, IINDERR, IINDWLC, IINFO, IIU, IM,
     $                   INDD, INDE, INDERR, INDILU, INDRW, INDWLC,
     $                   INDWORK, IPIL, IPIU, IPROC, ITMP, JTMP, LASTCL,
     $                   LENGTHI, LENGTHI2, LIWMIN, LWMIN, LWOPT,
     $                   MAXCLS, MQ00, MXCLSZ, MXDPTH, MYCOL, MYIL,
     $                   MYIU, MYPROC, MYROW, MZ, NB, NDEPTH, NEEDIL,
     $                   NEEDIU, NNP, NP00, NPCOL, NPROCS, NPROW,
     $                   NSPLIT, PARITY, RLENGTHI, RLENGTHI2, RSTARTI,
     $                   SIZE1, SIZE2, SRCCOL, SRCROW, STARTI, ZOFFSET

      DOUBLE PRECISION            PIVMIN, SAFMIN, SCALE, VLL, VUU, WL,
     $                            WU
*
*     .. Local Arrays ..
      INTEGER            IDUM1( 4 ), IDUM2( 4 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, INDXG2P, NUMROC, PJLAENV
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL            ICEIL, INDXG2P, LSAME, NUMROC, PDLAMCH,
     $                    PJLAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL            BLACS_GRIDINFO, CHK1MAT, DCOPY, DGEBR2D,
     $                    DGEBS2D, DGERV2D, DGESD2D, DLARRC_CV, DLASRT2,
     $                    DSTEGR2A_CV, DSTEGR2B_CV, DSTEGR2_CV, IGAMX2D,
     $                    IGEBR2D, IGEBS2D, IGERV2D, IGESD2D, IGSUM2D,
     $                    PCHK2MAT, PDLAEVSWP, PDSYNTRD, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, ICHAR, INT, MAX, MIN, MOD, SQRT
*     ..
*     .. Executable Statements ..
*

      INFO = 0

      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )

***********************************************************************
*
*     GET MACHINE PARAMETERS
*
***********************************************************************
      ICTXT = DESCZ( CTXT_ )
      SAFMIN = PDLAMCH( ICTXT, 'Safe minimum' )

***********************************************************************
*
*     Set up pointers into the WORK array
*     
***********************************************************************
*     Communication buffers for tridiagonal data, eigenvalues, errors 
      INDD = 1
      INDE = INDD + N+1
*     Workspace
      INDWORK = INDE + N+1

***********************************************************************
*
*     BLACS PROCESSOR GRID SETUP
*
***********************************************************************
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )


      NPROCS = NPROW * NPCOL
      MYPROC = MYROW * NPCOL + MYCOL
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 800+CTXT_ )
      ELSE IF( WANTZ ) THEN
         IF( ICTXT.NE.DESCZ( CTXT_ ) ) THEN
            INFO = -( 2100+CTXT_ )
         END IF
      END IF

***********************************************************************
*
*     COMPUTE REAL WORKSPACE
*
***********************************************************************
      IF ( ALLEIG ) THEN
         MZ = N
      ELSE IF ( INDEIG ) THEN
         MZ = IU - IL + 1
      ELSE
*        Take upper bound for VALEIG case
         MZ = N
      END IF
*     
      IF ( WANTZ ) THEN
*        Need to make sure enough space for PDLAEVSWP
         NB =  DESCZ( 6 )
         NP00 = NUMROC( N, NB, 0, 0, NPROW )
         MQ00 = NUMROC( MZ, NB, 0, 0, NPCOL )
         INDRW = INDWORK + MAX(18*N, NP00*MQ00)
         LWMIN = INDRW - 1 + (ICEIL(MZ, NPROCS) + 2)*N
      ELSE
         INDRW = INDWORK + 12*N
         LWMIN = INDRW - 1
      END IF
*     The code that validates the input requires 3 workspace entries
      LWMIN = MAX(3, LWMIN)
      LWOPT = LWMIN
*
      SIZE1 = INDRW - INDWORK

***********************************************************************
*
*     COMPUTE INTEGER WORKSPACE
*
***********************************************************************
      NNP = MAX( N, NPROCS+1, 4 )
      IF ( WANTZ ) THEN
        LIWMIN = 12*NNP + 2*N 
      ELSE
        LIWMIN = 10*NNP + 2*N
      END IF

***********************************************************************
*
*     Set up pointers into the IWORK array
*     
***********************************************************************
*     Pointer to eigenpair distribution over processors
      INDILU = LIWMIN - 2*NPROCS + 1            
      SIZE2 = INDILU - 2*N 

***********************************************************************
*
*     Test the input arguments.
*
***********************************************************************
      IF( INFO.EQ.0 ) THEN
         IF( WANTZ )
     $      CALL CHK1MAT( N, 4, N, 4, IZ, JZ, DESCZ, 21, INFO )
*
         IF( INFO.EQ.0 ) THEN
            IF( .NOT. WANTZ ) THEN
               INFO = -1
            ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
               INFO = -2
            ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
               INFO = -3
            ELSE IF( VALEIG .AND. N.GT.0 .AND. VU.LE.VL ) THEN
               INFO = -10
            ELSE IF( INDEIG .AND. ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) )
     $                THEN
               INFO = -11
            ELSE IF( INDEIG .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) )
     $                THEN
               INFO = -12
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -21
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -23
            END IF
         END IF
         IDUM2( 1 ) = 1
         IF( LOWER ) THEN
            IDUM1( 2 ) = ICHAR( 'L' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'U' )
         END IF
         IDUM2( 2 ) = 2
         IF( ALLEIG ) THEN
            IDUM1( 3 ) = ICHAR( 'A' )
         ELSE IF( INDEIG ) THEN
            IDUM1( 3 ) = ICHAR( 'I' )
         ELSE
            IDUM1( 3 ) = ICHAR( 'V' )
         END IF
         IDUM2( 3 ) = 3
         IF( LQUERY ) THEN
            IDUM1( 4 ) = -1
         ELSE
            IDUM1( 4 ) = 1
         END IF
         IDUM2( 4 ) = 4
         IF( WANTZ ) THEN
            IDUM1( 1 ) = ICHAR( 'V' )
            CALL PCHK2MAT( N, 4, N, 4, IZ, JZ, DESCZ, 8, N, 4, N, 4, IZ,
     $                     JZ, DESCZ, 21, 4, IDUM1, IDUM2, INFO )
         END IF
         WORK( 1 ) = DBLE( LWOPT )
         IWORK( 1 ) = LIWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDSTEGR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

***********************************************************************
*
*     Quick return if possible
*
***********************************************************************
      IF( N.EQ.0 ) THEN
         IF( WANTZ ) THEN
            NZ = 0
         END IF
         M = 0
         WORK( 1 ) = DBLE( LWOPT )
         IWORK( 1 ) = LIWMIN
         RETURN
      END IF

      IF( VALEIG ) THEN
         VLL = VL
         VUU = VU
      ELSE
         VLL = ZERO
         VUU = ZERO
      END IF

***********************************************************************
*
*     SET IIL, IIU
*
***********************************************************************
      IF ( ALLEIG ) THEN 
         IIL = 1
         IIU = N
      ELSE IF ( INDEIG ) THEN
         IIL = IL
         IIU = IU
      ELSE IF ( VALEIG ) THEN
         CALL DLARRC_CV('T', N, VLL, VUU, D, E,
     $                   SAFMIN, EIGCNT, IIL, IIU, INFO)
*        Refine upper bound N that was taken 
         MZ = EIGCNT
         IIL = IIL + 1
      ENDIF

      IF(MZ.EQ.0) THEN
         M = 0
         IF( WANTZ ) THEN
            NZ = 0
         END IF
         WORK( 1 ) = DBLE( LWOPT )
         IWORK( 1 ) = LIWMIN
         RETURN
      END IF

      MYIL = 0
      MYIU = 0
      M = 0
      IM = 0
      NDEPTH = 0
      MAXCLS = 1    

***********************************************************************
*
*     COMPUTE WORK ASSIGNMENTS
*
***********************************************************************

C     Each processor computes all work assignments
      CALL CMPIM2( IIL, IIU, NPROCS,
     $             IWORK(INDILU), IWORK(INDILU+NPROCS) )
C     find local work assignment
      MYIL = IWORK(INDILU+MYPROC)
      MYIU = IWORK(INDILU+NPROCS+MYPROC)


      ZOFFSET = MAX(0, MYIL - IIL - 1)
      FIRST = ( MYIL .EQ. IIL )


***********************************************************************
*
*     CALLS TO MRRR KERNEL
*
***********************************************************************
      IF(.NOT.WANTZ) THEN
         IINFO = 0
         IF ( MYIL.GT.0 ) THEN
            DOL = 1
            DOU = MYIU - MYIL + 1
            CALL DSTEGR2_CV( JOBZ, 'I', N, D,
     $                  E, VLL, VUU, MYIL, MYIU,
     $                  IM, W( 1 ), WORK( INDRW ), N, 
     $                  MYIU - MYIL + 1,
     $                  IWORK( 1 ), WORK( INDWORK ), SIZE1, 
     $                  IWORK( 2*N+1 ), SIZE2, 
     $                  DOL, DOU, ZOFFSET, IINFO )
*           DSTEGR2 zeroes out the entire W array, so we can't just give
*           it the part of W we need.  So here we copy the W entries into
*           their correct location
            DO 49 I = 1, IM
              W( MYIL-IIL+I ) = W( I )
 49         CONTINUE
*           W( MYIL ) is at W( MYIL - IIL + 1 )
*           W( X ) is at W(X - IIL + 1 )
         END IF
         IF (IINFO .NE. 0) THEN
            CALL PXERBLA( ICTXT, 'DSTEGR2', -IINFO )
            RETURN
         END IF
      ELSEIF ( WANTZ .AND. NPROCS.EQ.1 ) THEN
         IINFO = 0
         IF ( MYIL.GT.0 ) THEN
            DOL = MYIL - IIL + 1
            DOU = MYIU - IIL + 1
            CALL DSTEGR2_CV( JOBZ, 'I', N, D,
     $                  E, VLL, VUU, IIL, IIU,
     $                  IM, W( 1 ), WORK( INDRW ), N, 
     $                  N,
     $                  IWORK( 1 ), WORK( INDWORK ), SIZE1, 
     $                  IWORK( 2*N+1 ), SIZE2, DOL, DOU,
     $                  ZOFFSET, IINFO )
         ENDIF
         IF (IINFO .NE. 0) THEN
            CALL PXERBLA( ICTXT, 'DSTEGR2', -IINFO )
            RETURN
         END IF
      ELSEIF ( WANTZ ) THEN
*        Compute representations in parallel.
*        Share eigenvalue computation for root between all processors
*        Then compute the eigenvectors. 
         IINFO = 0
*        Part 1. compute root representations and root eigenvalues
         IF ( MYIL.GT.0 ) THEN
            DOL = MYIL - IIL + 1
            DOU = MYIU - IIL + 1
            CALL DSTEGR2A_CV( JOBZ, 'I', N, D,
     $                  E, VLL, VUU, IIL, IIU,
     $                  IM, W( 1 ), WORK( INDRW ), N, 
     $                  N, WORK( INDWORK ), SIZE1, 
     $                  IWORK( 2*N+1 ), SIZE2, DOL, 
     $                  DOU, NEEDIL, NEEDIU,
     $                  INDERR, NSPLIT, PIVMIN, SCALE, WL, WU,
     $                  IINFO )
         ENDIF
         IF (IINFO .NE. 0) THEN
            CALL PXERBLA( ICTXT, 'DSTEGR2A', -IINFO )
            RETURN
         END IF
*
         VSTART = .TRUE.
         FINISH = (MYIL.LE.0)
C        Part 2. Share eigenvalues and uncertainties between all processors
         IINDERR = INDWORK + INDERR - 1
*
*
         DOBCST = .TRUE.
         DOBCST = .FALSE.
         IF(DOBCST) THEN
*           First gather everything on the first processor.
*           Then use BROADCAST-based communication 
            DO 45 I = 2, NPROCS
               IF (MYPROC .EQ. (I - 1)) THEN
                  DSTROW = 0
                  DSTCOL = 0
                  STARTI = DOL
                  IWORK(1) = STARTI
                  IF(MYIL.GT.0) THEN
                     LENGTHI = MYIU - MYIL + 1
                  ELSE
                     LENGTHI = 0
                  ENDIF
                  IWORK(2) = LENGTHI
                  CALL IGESD2D( ICTXT, 2, 1, IWORK, 2, 
     $                    DSTROW, DSTCOL )
                  IF (( STARTI.GE.1 ) .AND. ( LENGTHI.GE.1 )) THEN
                     LENGTHI2 = 2*LENGTHI
*                    Copy eigenvalues into communication buffer
                     CALL DCOPY(LENGTHI,W( STARTI ),1,
     $                          WORK( INDD ), 1)                    
*                    Copy uncertainties into communication buffer
                     CALL DCOPY(LENGTHI,WORK( IINDERR+STARTI-1 ),1,
     $                          WORK( INDD+LENGTHI ), 1)                    
*                    send buffer
                     CALL DGESD2D( ICTXT, LENGTHI2, 
     $                    1, WORK( INDD ), LENGTHI2,
     $                    DSTROW, DSTCOL )
                  END IF
               ELSE IF (MYPROC .EQ. 0) THEN
                  SRCROW = (I-1) / NPCOL
                  SRCCOL = MOD(I-1, NPCOL)
                  CALL IGERV2D( ICTXT, 2, 1, IWORK, 2, 
     $                    SRCROW, SRCCOL )
                  STARTI = IWORK(1)
                  LENGTHI = IWORK(2)
                  IF (( STARTI.GE.1 ) .AND. ( LENGTHI.GE.1 )) THEN
                     LENGTHI2 = 2*LENGTHI
*                    receive buffer
                     CALL DGERV2D( ICTXT, LENGTHI2, 1,
     $                 WORK(INDD), LENGTHI2, SRCROW, SRCCOL )
*                    copy eigenvalues from communication buffer
                     CALL DCOPY( LENGTHI, WORK(INDD), 1,
     $                          W( STARTI ), 1)                    
*                    copy uncertainties (errors) from communication buffer
                     CALL DCOPY(LENGTHI,WORK(INDD+LENGTHI),1,
     $                          WORK( IINDERR+STARTI-1 ), 1)     
                  END IF
               END IF
  45        CONTINUE
            LENGTHI = IIU - IIL + 1
            LENGTHI2 = LENGTHI * 2
            IF (MYPROC .EQ. 0) THEN
*              Broadcast eigenvalues and errors to all processors
               CALL DCOPY(LENGTHI,W ,1, WORK( INDD ), 1)                 
               CALL DCOPY(LENGTHI,WORK( IINDERR ),1,
     $                          WORK( INDD+LENGTHI ), 1)                    
               CALL DGEBS2D( ICTXT, 'A', ' ', LENGTHI2, 1, 
     $              WORK(INDD), LENGTHI2 )
            ELSE
               SRCROW = 0
               SRCCOL = 0
               CALL DGEBR2D( ICTXT, 'A', ' ', LENGTHI2, 1,
     $             WORK(INDD), LENGTHI2, SRCROW, SRCCOL )
               CALL DCOPY( LENGTHI, WORK(INDD), 1, W, 1)
               CALL DCOPY(LENGTHI,WORK(INDD+LENGTHI),1,
     $                          WORK( IINDERR ), 1)                   
            END IF
         ELSE
*           Enable point2point communication between collaborators

*           Find collaborators of MYPROC            
            IF( (NPROCS.GT.1).AND.(MYIL.GT.0) ) THEN
               CALL CMPCOL( MYPROC, NPROCS, IIL, NEEDIL, NEEDIU, 
     $                   IWORK(INDILU), IWORK(INDILU+NPROCS),
     $                   COLBRT, FRSTCL, LASTCL )
            ELSE
               COLBRT = .FALSE.
            ENDIF

            IF(COLBRT) THEN
*              If the processor collaborates with others,
*              communicate information. 
               DO 47 IPROC = FRSTCL, LASTCL
                  IF (MYPROC .EQ. IPROC) THEN
                     STARTI = DOL
                     IWORK(1) = STARTI
                     LENGTHI = MYIU - MYIL + 1
                     IWORK(2) = LENGTHI
                     
                     IF ((STARTI.GE.1) .AND. (LENGTHI.GE.1)) THEN
*                       Copy eigenvalues into communication buffer
                        CALL DCOPY(LENGTHI,W( STARTI ),1,
     $                              WORK(INDD), 1)                    
*                       Copy uncertainties into communication buffer
                        CALL DCOPY(LENGTHI,
     $                          WORK( IINDERR+STARTI-1 ),1,
     $                          WORK(INDD+LENGTHI), 1)                    
                     ENDIF
                     DO 46 I = FRSTCL, LASTCL                      
                        IF(I.EQ.MYPROC) GOTO 46
                        DSTROW = I/ NPCOL
                        DSTCOL = MOD(I, NPCOL)
                        CALL IGESD2D( ICTXT, 2, 1, IWORK, 2, 
     $                             DSTROW, DSTCOL )
                        IF ((STARTI.GE.1) .AND. (LENGTHI.GE.1)) THEN
                           LENGTHI2 = 2*LENGTHI
*                          send buffer
                           CALL DGESD2D( ICTXT, LENGTHI2, 
     $                          1, WORK(INDD), LENGTHI2,
     $                          DSTROW, DSTCOL )
                        END IF
  46                 CONTINUE
                  ELSE
                     SRCROW = IPROC / NPCOL
                     SRCCOL = MOD(IPROC, NPCOL)
                     CALL IGERV2D( ICTXT, 2, 1, IWORK, 2, 
     $                             SRCROW, SRCCOL )
                     RSTARTI = IWORK(1)
                     RLENGTHI = IWORK(2)
                     IF ((RSTARTI.GE.1 ) .AND. (RLENGTHI.GE.1 )) THEN
                        RLENGTHI2 = 2*RLENGTHI
                        CALL DGERV2D( ICTXT, RLENGTHI2, 1,
     $                      WORK(INDE), RLENGTHI2,
     $                      SRCROW, SRCCOL )
*                       copy eigenvalues from communication buffer
                        CALL DCOPY( RLENGTHI, WORK(INDE), 1,
     $                          W( RSTARTI ), 1)                    
*                       copy uncertainties (errors) from communication buffer
                        CALL DCOPY(RLENGTHI,WORK(INDE+RLENGTHI),1,
     $                          WORK( IINDERR+RSTARTI-1 ), 1)                    
                     END IF
                  END IF
  47           CONTINUE
            ENDIF
         ENDIF

 100     CONTINUE
*
*        Part 3. Compute representation tree and eigenvectors
         IF ( MYIL.GT.0 ) THEN
            CALL DSTEGR2B_CV( JOBZ, N, D, E,
     $                  IM, W( 1 ), WORK( INDRW ), N, N,
     $                  IWORK( 1 ), WORK( INDWORK ), SIZE1, 
     $                  IWORK( 2*N+1 ), SIZE2, DOL, 
     $                  DOU, NEEDIL, NEEDIU, INDWLC,
     $                  PIVMIN, SCALE, WL, WU,
     $                  VSTART, FINISH, 
     $                  MAXCLS, NDEPTH, PARITY, ZOFFSET, IINFO )
            IINDWLC = INDWORK + INDWLC - 1
            IF(.NOT.FINISH) THEN
               IF((NEEDIL.LT.DOL).OR.(NEEDIU.GT.DOU)) THEN
                  CALL CMPCOL( MYPROC, NPROCS, IIL, NEEDIL, NEEDIU,
     $                 IWORK(INDILU), IWORK(INDILU+NPROCS),
     $                   COLBRT, FRSTCL, LASTCL )
               ELSE
                  COLBRT = .FALSE.
                  FRSTCL = MYPROC
                  LASTCL = MYPROC
               ENDIF

               IF(COLBRT) THEN
                  DO 147 IPROC = FRSTCL, LASTCL
                     IF (MYPROC .EQ. IPROC) THEN
                        STARTI = DOL
                        IWORK(1) = STARTI
                        IF(MYIL.GT.0) THEN
                           LENGTHI = MYIU - MYIL + 1
                        ELSE
                           LENGTHI = 0
                        ENDIF
                        IWORK(2) = LENGTHI
                        IF ((STARTI.GE.1).AND.(LENGTHI.GE.1)) THEN
*                          Copy eigenvalues into communication buffer
                           CALL DCOPY(LENGTHI,
     $                          WORK( IINDWLC+STARTI-1 ),1,
     $                          WORK(INDD), 1)                    
*                          Copy uncertainties into communication buffer
                           CALL DCOPY(LENGTHI,
     $                          WORK( IINDERR+STARTI-1 ),1,
     $                          WORK(INDD+LENGTHI), 1)                    
                        ENDIF
                     
                        DO 146 I = FRSTCL, LASTCL                      
                           IF(I.EQ.MYPROC) GOTO 146
                           DSTROW = I/ NPCOL
                           DSTCOL = MOD(I, NPCOL)
                           CALL IGESD2D( ICTXT, 2, 1, IWORK, 2, 
     $                             DSTROW, DSTCOL )
                           IF ((STARTI.GE.1).AND.(LENGTHI.GE.1)) THEN
                              LENGTHI2 = 2*LENGTHI

*                             send buffer
                              CALL DGESD2D( ICTXT, LENGTHI2, 
     $                             1, WORK(INDD), LENGTHI2,
     $                             DSTROW, DSTCOL )
                           END IF
 146                    CONTINUE
                     ELSE
                        SRCROW = IPROC / NPCOL
                        SRCCOL = MOD(IPROC, NPCOL)
                        CALL IGERV2D( ICTXT, 2, 1, IWORK, 2, 
     $                             SRCROW, SRCCOL )
                        RSTARTI = IWORK(1)
                        RLENGTHI = IWORK(2)
                        IF ((RSTARTI.GE.1).AND.(RLENGTHI.GE.1)) THEN
                           RLENGTHI2 = 2*RLENGTHI
                           CALL DGERV2D( ICTXT,RLENGTHI2, 1,
     $                         WORK(INDE),RLENGTHI2,
     $                         SRCROW, SRCCOL )
*                          copy eigenvalues from communication buffer
                           CALL DCOPY(RLENGTHI, WORK(INDE), 1,
     $                          WORK( IINDWLC+RSTARTI-1 ), 1)        
*                          copy uncertainties (errors) from communication buffer
                           CALL DCOPY(RLENGTHI,WORK(INDE+RLENGTHI),1,
     $                          WORK( IINDERR+RSTARTI-1 ), 1)            
                        END IF
                     END IF
 147              CONTINUE
               ENDIF
               GOTO 100         
            ENDIF
         ENDIF
         IF (IINFO .NE. 0) THEN
            CALL PXERBLA( ICTXT, 'DSTEGR2B', -IINFO )
            RETURN
         END IF
*
      ENDIF
*
***********************************************************************
*
*     MAIN PART ENDS HERE
*
***********************************************************************
*
***********************************************************************
*
*     ALLGATHER: EACH PROCESSOR SENDS ITS EIGENVALUES TO THE FIRST ONE,
*                THEN THE FIRST PROCESSOR BROADCASTS ALL EIGENVALUES
*
***********************************************************************
*
      DO 50 I = 2, NPROCS
         IF (MYPROC .EQ. (I - 1)) THEN
            DSTROW = 0
            DSTCOL = 0
            STARTI = MYIL - IIL + 1
            IWORK(1) = STARTI
            IF(MYIL.GT.0) THEN
               LENGTHI = MYIU - MYIL + 1
            ELSE
               LENGTHI = 0
            ENDIF
            IWORK(2) = LENGTHI
            CALL IGESD2D( ICTXT, 2, 1, IWORK, 2, 
     $                    DSTROW, DSTCOL )
            IF ((STARTI.GE.1).AND.(LENGTHI.GE.1)) THEN
               CALL DGESD2D( ICTXT, LENGTHI, 
     $              1, W( STARTI ), LENGTHI,
     $              DSTROW, DSTCOL )
            ENDIF
         ELSE IF (MYPROC .EQ. 0) THEN
            SRCROW = (I-1) / NPCOL
            SRCCOL = MOD(I-1, NPCOL)
            CALL IGERV2D( ICTXT, 2, 1, IWORK, 2, 
     $                    SRCROW, SRCCOL )
            STARTI = IWORK(1)
            LENGTHI = IWORK(2)
            IF ((STARTI.GE.1).AND.(LENGTHI.GE.1)) THEN
               CALL DGERV2D( ICTXT, LENGTHI, 1,
     $                 W( STARTI ), LENGTHI, SRCROW, SRCCOL )
            ENDIF
         ENDIF
   50 CONTINUE


*     Find maximum depth of tree among all processors
      MXDPTH = NDEPTH      
      CALL IGAMX2D( ICTXT, 'A', ' ', 1, 1, MXDPTH, 1, 
     $              ITMP, JTMP, -1, -1, -1 )
*     Find size of largest cluster
      MXCLSZ = MAXCLS
      CALL IGAMX2D( ICTXT, 'A', ' ', 1, 1, MXCLSZ, 1, 
     $              ITMP, JTMP, -1, -1, -1 )


*     Accumulate M from all processors
      M = IM
      CALL IGSUM2D( ICTXT, 'A', ' ', 1, 1, M, 1, -1, -1 )

*     Broadcast eigenvalues to all processors
      IF (MYPROC .EQ. 0) THEN
*        Send eigenvalues
         CALL DGEBS2D( ICTXT, 'A', ' ', M, 1, W, M )
      ELSE
         SRCROW = 0
         SRCCOL = 0
         CALL DGEBR2D( ICTXT, 'A', ' ', M, 1,
     $           W, M, SRCROW, SRCCOL )
      END IF
*
*     Sort the eigenvalues and keep permutation in IWORK to
*     sort the eigenvectors accordingly
*
      DO 160 I = 1, M
         IWORK( NPROCS+1+I ) = I
  160 CONTINUE
      CALL DLASRT2( 'I', M, W, IWORK( NPROCS+2 ), IINFO )
      IF (IINFO.NE.0) THEN
         CALL PXERBLA( ICTXT, 'DLASRT2', -IINFO )
         RETURN
      END IF

***********************************************************************
*
*     TRANSFORM Z FROM 1D WORKSPACE INTO 2D BLOCKCYCLIC STORAGE     
*
***********************************************************************
      IF ( WANTZ ) THEN
         DO 170 I = 1, M
            IWORK( M+NPROCS+1+IWORK( NPROCS+1+I ) ) = I
  170    CONTINUE
*        Store NVS in IWORK(1:NPROCS+1) for PDLAEVSWP
         IWORK( 1 ) = 0
         DO 180 I = 1, NPROCS
*           Find IL and IU for processor i-1
*           Has already been computed by CMPIM2 and stored
            IPIL = IWORK(INDILU+I-1)
            IPIU = IWORK(INDILU+NPROCS+I-1)
            IF (IPIL .EQ. 0) THEN
               IWORK( I + 1 ) = IWORK( I )
            ELSE
               IWORK( I + 1 ) = IWORK( I ) + IPIU - IPIL + 1
            ENDIF
  180    CONTINUE

         IF ( FIRST ) THEN
            CALL PDLAEVSWP(N, WORK( INDRW ), N, Z, IZ, JZ, 
     $       DESCZ, IWORK( 1 ), IWORK( NPROCS+M+2 ), WORK( INDWORK ), 
     $       INDRW - INDWORK )
         ELSE
            CALL PDLAEVSWP(N, WORK( INDRW + N ), N, Z, IZ, JZ, 
     $       DESCZ, IWORK( 1 ), IWORK( NPROCS+M+2 ), WORK( INDWORK ), 
     $       INDRW - INDWORK )
         END IF
*
         NZ = M
*
      END IF
*
      WORK( 1 ) = DBLE( LWOPT )
      IWORK( 1 ) = LIWMIN

      RETURN
*
*     End of PDSTEGR
*
      END
