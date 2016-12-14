!#begindoc
 
#ifndef WITH_UNION

#define diamatrix  anymatrix

#endif

MODULE m_invbdia

CONTAINS
 
SUBROUTINE invbdia (n, Ad)

USE m_dump
USE m_dgedi
USE m_dgeco
USE m_build

INTEGER						, INTENT(IN)            :: n
TYPE (diamatrix)				, POINTER         	:: Ad

!     Compute the Inverse of a Block-Diagonal, stored by column,
!     matrix Ad and store the result in Ad:
!        Ad := inv(Ad)

!     Arguments:
!     ==========
!     N        	i   Number of rows/columns in the block-diagonal
!                   matrix Ad.
!     Ad%blksiz	i   Block Size: number of rows/columns per diagonal block.
!                  'N' should be an integer multiple of 'Ad%blksiz'!
!     Ad%com     io  In:  Ad_in, elements of the diagonal blocks in
!                       matrix Ad.
!                   Out: elements of the diagonal blocks of the
!                       matrix  inv(Ad_in).

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'invbdia'

!     Local Variables:
!     ================
!     NBlk   Number of blocks in block-diagonal matrix  Ad.

INTEGER				:: ier
INTEGER 			:: nblk, firrow
INTEGER, DIMENSION(1:Ad%blksiz)	:: temp
DOUBLE PRECISION 		:: rcond
DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Adco

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Compute and check the value of  NBlk:
nblk = n / Ad%blksiz

IF (nblk*Ad%blksiz /= n) THEN
!        Error if 'N' and 'NBlk' do not satisfy:
!        "N is an integer multiple of Ad%blksiz"
  PRINT '( /, A, 2X, A, /, 2(3X, A, I10, /) )' ,  &
      'Internal error in', rounam, '!',  &
      'Block size:                           ', Ad%blksiz,  &
      'inconsistent with number rows/columns:', n
  CALL dump(__FILE__,__LINE__,'Internal error')
END IF

IF (Ad%blksiz /= 1) THEN
!        Size of diagonal blocks > 1:
  
!        For each block
  DO firrow = 1, nblk
!           Factor the matrix and estimate condition number:
    Adco => Ad%com( 1:Ad%blksiz, (firrow-1)*Ad%blksiz+1:firrow*Ad%blksiz )
    CALL dgeco (Adco, Ad%blksiz, temp, rcond)
    
!           Error, if matrix exactly singular or Rcond underflows
    IF (rcond <= 0.0D0) THEN
      PRINT '( /, A, 2X, A, /, 3X, A, / )' , 'Fatal error occurred in', rounam,  &
     'Singular block/element in (block-)diagonal matrix!'
      CALL dump(__FILE__,__LINE__,'Matrix exactly singular or Rcond underflows')
    END IF
    
!           Compute the inverse of a matrix

    CALL dgedi (Adco, Ad%blksiz, temp)
  END DO
  
ELSE
!        Size of diagonal blocks = 1:
  DO firrow = 1, nblk
    IF (DABS(Ad%com(1,firrow)) <= 0.0D0) THEN 
      PRINT '( /, A, 2X, A, /, 3X, A, / )' , 'Fatal error occurred in', rounam,  &
     'Singular block/element in (block-)diagonal matrix!'
      CALL dump(__FILE__,__LINE__,'Fatal error')
    END IF
    Ad%com(1,firrow) = 1.0D0/Ad%com(1,firrow)
  END DO
END IF

END SUBROUTINE invbdia

END MODULE