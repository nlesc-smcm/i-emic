!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

#ifdef MKL

#else

#define asum(x) SUM(ABS(x))

#endif

MODULE m_scalmat

CONTAINS
 
SUBROUTINE scalmat (a, scafctr)

USE m_glbpars
USE m_prcpars
USE m_vispars
USE m_wrtmt
USE m_ioerrmsg
USE m_build
USE m_dump

#ifdef MKL
USE mkl95_precision, ONLY: wp => dp
USE mkl95_blas, ONLY: asum
#endif

TYPE (csrmatrix)			, POINTER		:: a
DOUBLE PRECISION, DIMENSION(1:a%n)	, INTENT(OUT)           :: scafctr

!     Scales the A%N x A%N matrix A, stored in CSR format, according to
!     the value of 'ScaRow' in common block /prcpars/.
!     The scale factors are stored into 'ScaFctr', so that
!        A_out := diag(ScaFctr) A_in

!     Arguments:
!     ==========
!     A%N        i   Number of rows/columns in matrix A.
!     a%beg     i   a%beg(i): index in 'a%jco' and 'a%co' of the first
!                  nonzero element in row i of A.
!     a%jco     i   a%jco(nz): column number of the non-zero element
!                  a%co(nz).
!     a%co      io  a%co(nz): value of non-zero element of  A.
!                  In : non-zero value of element in matrix  A.
!                  Out: scaled value of this element.
!     ScaFctr  o   ScaFctr(i), 1 <= i <= A%N, the scaling factor applied
!                  to row 'i' of the matrix A.

!#enddoc

!     Global Parameters:
!     ==================

!     Local Variables:
!     ================
INTEGER 		:: i, nz, p, q, ier
DOUBLE PRECISION 	:: diaval, fact, MaxValA, minnz, rowsum, scafct

CHARACTER (LEN=*), PARAMETER :: rounam = 'scalmat'

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif


IF (visasc) THEN
!        Visualisation of scaled matrix  A  required:
  
  IF (outlev >= 1) THEN
    PRINT '(A, /)' , 'Visualize the scaled matrix with:  vsm APrescaled'
  END IF
  CALL  wrtmt ('APrescaled', .true., a, 2, ier)
  IF (ier /= 0) CALL ioerrmsg ('APrescaled', 'wrtmt' , rounam, ier, 0)
END IF


!     Determine matrix element with maximum absolute value:
!MaxValA = ABS(a%co(idamax (csrnnz(a), a%co(1:csrnnz(a)) )))
MaxValA = MAXVAL( ABS( a%co(1:csrnnz(a)) ) )
minnz  = neglgbl * MaxValA

IF (MaxValA <= 1.0D-96) THEN
  PRINT '(A, 2X, A, A, /, 3X, A, 1P, E12.5, A)',  &
      'Fatal error in:', rounam, '!',  &
      'Maximum absolute value of matrix elements,', MaxValA, ', is too small.'
  CALL dump(__FILE__,__LINE__, 'Fatal error in, maximum absolute value of matrix elements is too small!')
ELSE IF (MaxValA >= 1.0D+96) THEN
  PRINT '(A, 2X, A, A, /, 3X, A, 1P, E12.5, A)',  &
      'Fatal error in:', rounam, '!',  &
      'Maximum absolute value of matrix elements,', MaxValA, ', is too large.'
  CALL dump(__FILE__,__LINE__, 'Fatal error in, maximum absolute value of matrix elements is too large!')
END IF

IF (.NOT. scarow) THEN

  ScaFct = 2**( FLOOR( LOG(MaxValA)/LOG(2.0D0 ) ) )

! {  0 < ScaFct <= MaxValA < 2*ScaFct  }

  scafct = 1.0D0 / scafct

! {  1 <= ScaFct*MaxValA < 2  }

END IF

DO i = 1, a%n

  p = a%beg(i)
  q = a%beg(i+1)-1

! Determine value of diagonal element:

  diaval = minnz
  DO nz = p, q
    IF (a%jco(nz) == i) THEN
      diaval = a%co(nz)
      EXIT
    END IF
  END DO
  
  IF (scarow) THEN

!   Compute  ScaFct = 1 / (MinNz + SUM(j :: ABS(A(i,j)))):

    scafct = 1.0D0 / ( minnz + asum( a%co(p:q) ) )
  END IF
  
!  Modify scale factor to ensure positive diagonal:

  IF (diaval < 0) THEN
    fact = -scafct
  ELSE
    fact = scafct
  END IF
  scafctr(i) = fact
  
! Apply the row scaling:
  
  a%co(p:q) = fact * a%co(p:q)

END DO


IF (visasc) THEN
!        Visualisation of scaled matrix  A  required:
  
  IF (outlev >= 1) THEN
    PRINT '(A, /)' , 'Visualize the scaled matrix with:  vsm Ascaled'
  END IF
  CALL  wrtmt ('Ascaled', .true., a, 2, ier)
  IF (ier /= 0) CALL ioerrmsg ('Ascaled', 'wrtbmt' , rounam, ier, 0)
END IF


END SUBROUTINE scalmat

END MODULE