!#begindoc
 
#ifndef WITH_UNION

#define scdematrix  anymatrix

#endif

MODULE m_visscd

CONTAINS

SUBROUTINE visscd (message, filenm, xs)

USE m_build
USE m_glbpars
USE m_wrtmtd
USE m_ioerrmsg

CHARACTER (LEN=*)		, INTENT(IN)		:: message
CHARACTER (LEN=*)		, INTENT(IN)            :: filenm
TYPE (scdematrix)		, POINTER		:: xs

!     Visualize the SCD type matrix [xs].

!     Writes the SCD type matrix [xs] to the binary file 'filenm'
!     in CSR format and writes the text
!        'message'  vsm 'filenm
!     to the standard output file.

!     Arguments:
!     ==========
!     message  i   The message string
!     filenm   i   File name string
!     xs      i   Location of descriptor for an SCD type
!                  matrix, containing the Schur-complement of
!                  A_11  in  A.

!#enddoc

!     Local Parameters:
!     =================
CHARACTER (LEN=*), PARAMETER :: rounam = 'visscd'

  IF (outlev >= 1) THEN
    PRINT '(A, 2X, A, X, A, /)' ,  message, 'vsm', filenm
  END IF
  
  CALL wrtmtd (filenm, .true., xs%offd, xs%dia, 2)
	  
END SUBROUTINE visscd

END MODULE