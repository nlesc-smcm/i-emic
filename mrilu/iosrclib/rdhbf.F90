!#begindoc

#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define diamatrix  anymatrix

#endif
 
MODULE m_rdHBF

CONTAINS

SUBROUTINE rdhbf (filenm, binary,a,b,x0,x,rhs,guess,exact)

USE m_dump
USE m_getunit
USE m_fstrlen
USE m_build
USE m_getunit
USE m_wacsr

CHARACTER (LEN=*)			, INTENT(IN)            :: filenm
LOGICAL					, INTENT(IN)		:: binary
LOGICAL					, INTENT(out)		:: rhs,guess,exact
TYPE (csrmatrix)			, POINTER		:: a
DOUBLE PRECISION, DIMENSION(:)		, POINTER		:: b,x0,x

!     Reads a vector V from the ASCII file 'filenm'.
!     The representation of the vector in the file is
!        N, v
!     Numbers in the file should be separated by at least one space or
!     end-of-line!
!     The order of the numbers in the file are
!        1 * integer:                the length of vector  V
!        End-of-line
!        N * double precision:       all values of V(1:N)

!     Arguments:
!     ==========
!     filenm   i   Filename
!     N        o   The actual length of the vector V.
!     v        o   v(i), 1 <= i <= N, i-th element of vector V.
!     ier      o   The error code,
!                  =  0   No error
!                  = -3   Error during closing, opening or rewinding
!                         file 'filenm'.
!                  = -4   Unexpected End Of File encountered.
!                  = -5   Error during reading  'filenm'.
!                  =-12   No logical unit available to read from
!                         'filenm'.



!#enddoc

!     1997-11-12  ,  Doeke de Vries
!     2003-03-05  ,  Last update (Doeke de Vries)

!     ================================================================
!     ... SAMPLE CODE FOR READING A GENERAL SPARSE MATRIX, POSSIBLY
!         WITH RIGHT-HAND SIDE VECTORS
!     ================================================================

      CHARACTER      TITLE*72 , KEY*8    , MXTYPE*3 , RHSTYP*3
      CHARACTER      PTRFMT*16, INDFMT*16, VALFMT*20, RHSFMT*20

      INTEGER        I,ier,n,TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD
      INTEGER           NROW  , NCOL  , NNZERO, NELTVL
      INTEGER           NRHS  , NRHSIX, NRHSVL, NGUESS, NEXACT

!      INTEGER        POINTR (*), ROWIND (*), RHSPTR (*), RHSIND(*)
      INTEGER        RHSPTR (1), RHSIND(1)

!      REAL           VALUES (*) , RHSVAL (*), XEXACT (*), SGUESS (*)
      REAL           VALUES (1) , RHSVAL (1), XEXACT (1), SGUESS (1)
      INTEGER     LUNIT
!     Get an unopened logical unit number in LUNIT
      CALL getunit (LUNIT, ier)
      IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Could not find free unit')     
      OPEN (LUNIT, FILE = filenm, STATUS = 'OLD', ERR = 997)
 997  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Error opening file')     
!     ------------------------
!     ... READ IN HEADER BLOCK
!     ------------------------

      READ ( LUNIT, 1000 ) TITLE , KEY   ,&
                          TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD, &
                          MXTYPE, NROW  , NCOL  , NNZERO, NELTVL, &
                          PTRFMT, INDFMT, VALFMT, RHSFMT
      CALL wacsr (NCOL, NNZERO, a)
      
      IF  ( RHSCRD .GT. 0 )&
          READ ( LUNIT, 1001 ) RHSTYP, NRHS, NRHSIX
      IF (NRHS.GE.1) THEN
         PRINT *, 'NRHS =', NRHS
      ENDIF
      IF (NRHS.GE.1) THEN
        IF (RHSTYP(1:1) .NE. 'F' ) THEN 
          NRHS=1
        ENDIF
      ENDIF
      n=NCOL
      rhs=.FALSE.
      guess=.FALSE.
      exact=.FALSE.   

 1000 FORMAT ( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )
 1001 FORMAT ( A3, 11X, 2I14 )

!     -------------------------
!     ... READ MATRIX STRUCTURE
!     -------------------------

      READ ( LUNIT, PTRFMT ) ( a%beg(I), I = 1, NCOL+1 )

      READ ( LUNIT, INDFMT ) ( a%jco (I), I = 1, NNZERO )

      IF  ( VALCRD .GT. 0 )  THEN

!         ----------------------
!         ... READ MATRIX VALUES
!         ----------------------

          IF  ( MXTYPE (3:3) .EQ. 'A' )  THEN
              READ ( LUNIT, VALFMT ) ( a%co (I), I = 1, NNZERO )
          ELSE
              READ ( LUNIT, VALFMT ) ( a%co (I), I = 1, NELTVL )
          ENDIF

!         -------------------------
!         ... READ RIGHT-HAND SIDES
!         -------------------------

          IF  ( NRHS .GT. 0 )  THEN

              IF  ( RHSTYP(1:1) .EQ. 'F' ) THEN
                  rhs=.TRUE.
                  ALLOCATE( b(1:n), STAT=ier)
                  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
!                 -------------------------------
!                 ... READ DENSE RIGHT-HAND SIDES
!                 -------------------------------

                  NRHSVL = NROW * NRHS
                  READ ( LUNIT, RHSFMT ) ( b (I), I = 1, NRHSVL )

              ELSE

!                 ---------------------------------------------
!                 ... READ SPARSE OR ELEMENTAL RIGHT-HAND SIDES
!                 ---------------------------------------------


                  IF (MXTYPE(3:3) .EQ. 'A') THEN

!                    ------------------------------------------------
!                    ... SPARSE RIGHT-HAND SIDES - READ POINTER ARRAY
!1                    ------------------------------------------------

                     READ (LUNIT, PTRFMT) ( RHSPTR (I), I = 1, NRHS+1 )

!                    ----------------------------------------
!                    ... READ SPARSITY PATTERN FOR RIGHT-HAND
!                        SIDES
!                    ----------------------------------------

                     READ (LUNIT, INDFMT) ( RHSIND (I), I = 1, NRHSIX )

!                    --------------------------------------
!                    ... READ SPARSE RIGHT-HAND SIDE VALUES
!                    --------------------------------------

                     READ (LUNIT, RHSFMT) ( RHSVAL (I), I = 1, NRHSIX )

                  ELSE

!                    -----------------------------------
!                    ... READ ELEMENTAL RIGHT-HAND SIDES
!                    -----------------------------------

                     NRHSVL = NNZERO * NRHS
                     READ (LUNIT, RHSFMT) ( RHSVAL (I), I = 1, NRHSVL )

                  ENDIF

              END IF

              IF  ( RHSTYP(2:2) .EQ. 'G' ) THEN
              guess=.TRUE.
              ALLOCATE( x0(1:n), STAT=ier)
              IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
!                 -------------------------
!                 ... READ STARTING GUESSES
!                 -------------------------

                 NGUESS = NROW * NRHS
                 READ (LUNIT, RHSFMT) ( x0 (I), I = 1, NGUESS )

              END IF

              IF  ( RHSTYP(3:3) .EQ. 'X' ) THEN
                 exact=.TRUE.
!                 -------------------------
!                 ... READ SOLUTION VECTORS
!                 -------------------------
                 ALLOCATE( x(1:n), STAT=ier)
                 IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
                 NEXACT = NROW * NRHS
                 READ (LUNIT, RHSFMT) ( x(I), I = 1, NEXACT )
              END IF
          END IF
      END IF
END SUBROUTINE rdhbf
END MODULE