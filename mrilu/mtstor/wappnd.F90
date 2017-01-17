!#begindoc
 
#ifndef WITH_UNION

#define mlpmatrix  anymatrix
#define partmatrix  anymatrix

#endif

MODULE m_wappnd

CONTAINS

SUBROUTINE wappnd (Part, x)

USE m_build

TYPE (partmatrix)	, POINTER	:: Part
TYPE (mlpmatrix)	, POINTER       :: x

!     Append descriptor of partition 'Part' to the list of partitions
!     in the MLP matrix 'x'.

!     Arguments:
!     ==========
!     Part   	i   Location of the Partition descriptor.
!     x    	i   Location of descriptor of the Multi-Level
!                   Partitioned matrix P.  Also the head of the doubly
!                   linked list of partitions.

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'wappnd'

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam

PRINT '(A)', 'Allocating PART-matrix'
#endif

!     Insert the descriptor 'Part' at the end of the doubly linked:

IF ( ASSOCIATED(Part)) THEN

  IF ( ASSOCIATED(x%first) ) THEN
!   the second or later partition  
    Part%prev => x%last 
    Part%prev%next  => Part
  ELSE
!   the first partition
    NULLIFY(Part%prev)
    x%first => Part
  END IF
  NULLIFY(Part%next)
  x%last => Part

END IF

RETURN

END SUBROUTINE wappnd

END MODULE