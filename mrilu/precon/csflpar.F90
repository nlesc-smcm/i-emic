!#begindoc
 
#ifndef WITH_UNION

#define scdematrix  anymatrix
#define partmatrix  anymatrix

#endif

MODULE m_csflpar

CONTAINS

SUBROUTINE csflpar (NEqDon, NEqNotDon, S, Part)

USE m_dump
USE m_build
USE m_wapffp
USE m_wcompr
USE m_wfree
USE m_prcpars
USE m_scd2fm

#ifdef WITH_ATLAS 

EXTERNAL :: dgetrf

#else

USE m_dgetrf

#endif

INTEGER			, INTENT(IN)		:: NEqDon
INTEGER			, INTENT(IN)            :: NEqNotDon
TYPE (scdematrix)	, POINTER		:: S
TYPE (partmatrix)	, POINTER		:: Part

!     Compute and Store Full Last Partition into a new area .
!     The storage of the matrix  S  is released.

!     Arguments:
!     ==========
!     NEqDon   i   Number of rows/columns that have been entered in the
!                  Multilevel Partitioned LDU matrix.
!     NEqNotDon   i   Number of rows/columns in Schur-complement matrix
!                  that has to be factored.
!     S      io  In:  Location of the
!                       Schur-complement matrix  S.
!                  Out: -- undefined
!     Part   o   Location of descriptor of the newly
!                  allocated last partition with the L, D and U factors.

!#enddoc

!     Local variables:
!     ================
INTEGER 					:: i, inxmin, ier
DOUBLE PRECISION 				:: abselm, minelm, maxelm
DOUBLE PRECISION, DIMENSION(:,:), POINTER 	:: U

CHARACTER (LEN=*), PARAMETER :: rounam = 'csflpar'

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Get the constituent parts of the representation of matrix S:

!     Allocate descriptors and segments for a Partition with Full Factors and Pivots and append the
!     partition to the list in the MLP matrix  P [P]:
!     (Descriptor has been initialised!)
CALL wapffp (NEqNotDon, NEqDon, Part)

!     Copy the Schur-complement S, in SCD format, into the Full Matrix
!     of the PFFP partition:
CALL scd2fm (NEqNotDon, S, Part%fm)

!     Release workspace occupied by matrix S:
CALL scdefree(S)

!     Get the base index of the full matrix:
U => Part%fm%com

!     Compute the LU factorization of the matrix  B.  Store the factors
!     L  and  U  in  B  and the pivot indices in Piv:

#ifdef WITH_ATLAS 

CALL dgetrf(Part%N,Part%N,U,Part%N,Part%piv,ier)
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Failure in dgetref')

#else

CALL dgetrf (U, Part%piv, ier)
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Failure in dgetref')

#endif


IF (Part%N > 1) THEN
  
!        Check if upper triangular matrix U is nearly singular.
!        Compute:
!           maxelm = MAX(i : 1<=i<=Part%N : ABS(U(i,i)) )
!           minelm = MIN(i : 1<=i< Part%N : ABS(U(i,i)) )
  
  abselm = ABS(U(1,1))
  
  inxmin = 1
  minelm = abselm
  maxelm = abselm
  
  DO i = 2, Part%N - 1
    abselm = ABS(U(i,i))
    IF (abselm < minelm) THEN
      inxmin = i
      minelm = abselm
    ELSE IF (abselm > maxelm) THEN
      maxelm = abselm
    END IF
  END DO

  abselm = ABS(U(Part%N,Part%N))
  IF (abselm > maxelm)  maxelm = abselm
  
!        abselm = ABS(U(Part%N,Part%N))
!        minelm = MIN(i: 1<=i< Part%N: ABS(U(i,i))) = U(inxmin,inxmin)
!        maxelm = MAX(i: 1<=i<=Part%N: ABS(U(i,i)))
  
  IF (minelm <= lutol * maxelm) THEN
!           U(inxmin,inxmin) is (nearly) zero; 1 <= inxmin < Part%N.
!           Never allowed!
    PRINT '(/, A, X, A, /, 3X, A, /, 3X, I11, A, /)',  &
        'Error detected after  DGETRF  called in', rounam,  &
        'After LU-factorisation of last submatrix:',  &
        inxmin, '-th diagonal element of U is (nearly) zero!'
    CALL dump(__FILE__,__LINE__,'Failure')
  END IF
  
  IF (abselm <= lutol * maxelm) THEN
!           U(Part%N,Part%N) is nearly zero.
    IF (.NOT. singlu) THEN
      PRINT '(/, A, X, A, /, 3X, A, /, 3X, A, /)',  &
          'Error detected after  DGETRF  called in', rounam,  &
          'After LU-factorisation of last submatrix:',  &
          'Last diagonal element of U is (nearly) zero!'
      CALL dump(__FILE__,__LINE__,'Failure')
    END IF
    
!           (Nearly) Singular matrix U is allowed.
    PRINT '(/, A, 2(/, 3X, A), 1P, E12.5, /)',  &
        'Warning after LU-factorisation of last submatrix:',  &
        'Last diagonal element of factor U is (nearly) zero!',  &
        'Value of this element is changed into', maxelm
    
    U(Part%N,Part%N) = maxelm
  END IF
  
ENDIF

!     Compress the storage by shifting the
!     representation of 'Part' downwards and adjust the required size:
CALL wcompr (parttoany(Part))

!     End of  csflpar
END SUBROUTINE csflpar

END MODULE