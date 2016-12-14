!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

MODULE m_csrresize

CONTAINS

SUBROUTINE csrresize (x, maxnnz)

USE m_dump
USE m_build

TYPE (csrmatrix), POINTER               :: x
INTEGER, INTENT(IN)               	:: maxnnz

!     Resizes the storage occupied by the descriptor, indicated by 'x', with
!     all the referenced segments from this descriptor. The new size is maxnnz non
!     zeros.

!     Arguments:
!     ==========
!     x   	io  The matrix descriptor.
!     maxnnz    i   New number of non zeros

!#enddoc

!     Local Parameters:                                 :
!     =================

CHARACTER (LEN=*), PARAMETER :: rounam = 'csrresize'

!     Local Variables:
!     ================

INTEGER 				:: i, n, ier
INTEGER, DIMENSION(:), POINTER 		:: iarray
DOUBLE PRECISION, DIMENSION(:), POINTER :: darray, dold

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)', 'Entry:', rounam
#endif

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

!x%nnz = csrnnz(x)
!PRINT '(A, 3I, A)', 'Entry resize:', x%nnz , csrnnz(x), SUM( x%beg(2:x%n+1) - x%beg(1:x%n) )
IF (maxnnz < x%nnz) CALL dump(__FILE__,__LINE__,'Too small size')

IF (maxnnz /= UBOUND(x%co,1)) THEN    
  ALLOCATE( darray(1:maxnnz), STAT=ier)
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

! BUG: the next statement leeds in openmp mode for large array's to a segmentation violation.

! darray=x%co(1:x%nnz)

! alternative implementation

  DO i=1,x%nnz
    darray(i)=x%co(i)
  ENDDO

  DEALLOCATE( x%co, STAT=ier)
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
  x%co=>darray
END IF

IF (maxnnz /= UBOUND(x%jco,1)) THEN    
  ALLOCATE( iarray(1:maxnnz), STAT=ier)
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

! BUG: the next statement leeds in openmp mode for large array's to a segmentation violation.

! iarray=x%jco(1:x%nnz)

! alternative implementation

  DO i=1,x%nnz
    iarray(i)=x%jco(i)
  ENDDO

  DEALLOCATE( x%jco, STAT=ier)
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
  x%jco=>iarray
END IF

!     Normal Return:

RETURN
      
END SUBROUTINE csrresize

END MODULE
