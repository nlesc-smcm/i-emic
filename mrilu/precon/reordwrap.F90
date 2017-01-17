!#begindoc
 
#ifndef WITH_UNION

#define scdematrix  anymatrix
#define csrdmatrix  anymatrix

#endif

MODULE m_reordwrap

CONTAINS

SUBROUTINE reordwrap (n, Nupp, BlkSiz, S, B, Permut, CSpace, RSpace)

USE m_dump
USE m_build
USE m_wacsr
USE m_wacsrd
USE m_wcompr
USE m_glbpars
USE m_prcpars
USE m_vispars
USE m_xfnminr
USE m_fstrlen
USE m_wrtmtd
USE m_ioerrmsg
USE m_hernumsch
USE m_lumpspace
USE m_lumpdrop
USE m_order
USE m_permbdi
USE m_iperv
USE m_dperv

INTEGER					, INTENT(IN)		:: n
INTEGER					, INTENT(OUT)           :: Nupp
INTEGER					, INTENT(IN)            :: BlkSiz
TYPE (scdematrix)			, POINTER		:: S
TYPE (csrdmatrix)			, POINTER		:: B
INTEGER, DIMENSION(:) 			, POINTER	        :: Permut
DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: CSpace
DOUBLE PRECISION, DIMENSION(1:n)	, INTENT(IN OUT)        :: RSpace

!     Computes a new order of rows/columns in the  N x N  sub-
!     matrix  S.  Stores the reordered matrix into the newly created
!     matrix B, dropping or lumping elements.  The matrix  B  is stored.
!     The reordering is reflected in the update of the permutation
!     vector 'Permut', and the (used) lump spaces are updated in
!     'CSpace' and 'RSpace'.

!     Arguments:
!     ==========
!     N        	i   Number of rows/columns in matrix S and matrix B.
!     Nupp     	o   Number of rows in the upper partition of matrix B.
!     S      	i   Location of the descriptor of the Schur-complement matrix  S.
!     B      	o   Location of the reordered matrix  B, stored in CSRD format.
!     Permut   	io  Segment of the permutation vector relevant for the
!                   lower right part of the multilevel preconditioner.
!     CSpace   	io  The (updated) Column lump space from the 1st Schur-complement if  CLSOnce,  or
!                   The used Column lump space.
!     RSpace    io  The (updated) Row lump space from the 1st Schur-complement if  CLSOnce,  or
!                   The used Row lump space, if  .NOT. CLSOnce.

!#enddoc

!     Local Parameters:
!     =================

CHARACTER (LEN=*), PARAMETER :: rounam = 'reordwrap'

INTEGER, ALLOCATABLE, DIMENSION(:)		:: NewPer
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: NewCS, NewRS, ActCS, ActRS

!     Local Variables:
!     ================
!     ActCS  the actual Column lump space.
!     ActRS  the actual Row lump space.
!     NewCS  the new Column lump space, and
!              of the initial value of the actual Column lump space.
!     NewRS  the new Row lump space, and
!              of the initial value of the actual Row lump space.
!     NewPer of the New Permutation  NewPer.

CHARACTER (LEN=20) 				:: rofnm
INTEGER 					:: fnmlen, ier
INTEGER 					:: nblock

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

ALLOCATE( NewPer(1:n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( NewCS(1:n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( NewRS(1:n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( ActCS(1:n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( ActRS(1:n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')


nblock = n / BlkSiz


IF (.NOT. clsonce) THEN
  
!        Set up Row- and Column Lump space for Preconditioner:
  IF (S%n /= N) THEN
    PRINT '( /, A, I10, A, I10 /)' , 'Number rows/columns in block diagonal: ', S%n, ' differs from number rows/columns in matrix: ', N
    CALL dump(__FILE__,__LINE__,'Internal error')
  END IF
  
  CALL lumpspace (S, NewCS, NewRS)
  
!        Compute the actual Column- and Row- lump spaces:
  
  ActCS = MAX(NewCS - CSpace, NewCS * nlsfctr  )
  NewCS = ActCS
  ActRS = MAX(NewRS - RSpace, NewRS * nlsfctr  )
  NewRS = ActRS
END IF

!     Calculate a new reordering of the Schur-complement  S  and
!     store the corresponding permutation vector in  NewPer.
IF (clsonce) THEN
  CALL order (BlkSiz, Nupp, S%offd, NewPer, CSpace, RSpace )
ELSE
  CALL order (BlkSiz, Nupp, S%offd, NewPer, ActCS, ActRS )
END IF

!     Test for fatal error
IF (Nupp < 1) CALL dump(__FILE__,__LINE__,'Too few unknowns can be eliminated approximately. Increase drop tolerance.')

!     Store part of the reordered off-diagonal part of S, into the
!     off-diagonal part of B:
!        B := P A inv(P) - R,
!     drop the small elements of B, or lump those elements onto Sd,
!     and return location of boundary between lower partition blocks.

!     Lump or drop the "small" off-diagonal elements in the new left-
!     and new upper part.
!     Update the lump space for the rows/columns in the new left-
!     and new upper part, indicated by  NewPer(1:Nupp):
IF (clsonce) THEN
  CALL lumpdrop (Nupp, NewPer, S, CSpace, RSpace )
ELSE
!        { .NOT. CLSOnce }
  CALL lumpdrop (Nupp, NewPer, S, ActCS, ActRS )
  
!        Update the Used Column- and Row- lump spaces:
  CSpace = CSpace + (NewCS - ActCS)
  RSpace = RSpace + (NewRS - ActRS)
  
END IF


!     Update the global renumbering, in Permut(1:N),
!     of the original equations and reorder also the lump spaces,
!     in  CSpace(1:N)  and  RSpace(1:N).
!     All permutations are based on 'NewPer'.
CALL iperv (.false., NewPer, Permut, Permut)
CALL dperv (.false., NewPer, CSpace, CSpace)
CALL dperv (.false., NewPer, RSpace, RSpace)

!     Request workspace for a Schur Complement, (block-) Diagonal
!     Extracted, matrix B:
CALL wacsrd (n, BlkSiz, csrnnz(S%offd), B)

!     Place the reordered and pruned matrix into  B:
CALL hernumsch (Nupp, NewPer, S%offd, B )

!     Store reordered Block Diagonal part of S, Sd, with the lumped
!     elements into the (block-) Block Diagonal part of B, Bd:
!        Bd := P Sd inv(P)
CALL permbdi (n, NewPer, S%dia, B%dia )

DEALLOCATE( NewPer, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( NewCS, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( NewRS, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( ActCS, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( ActRS, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

IF (visnro > nrowrtn) THEN
!        Visualisation of reordered matrix required
  
  nrowrtn = nrowrtn + 1
  
  CALL xfnminr ('Breord', 2, nrowrtn, rofnm)
  fnmlen = fstrlen(rofnm)
  IF (outlev >= 1) THEN
    PRINT '(A, X, A)' , 'Visualize the reordered matrix with:  vsm', rofnm(1:fnmlen)
  END IF
  CALL wrtmtd (rofnm(1:fnmlen), .true., B%offd, B%dia, 2)
END IF


!     Compress the storage 'B' adjust the required size:
CALL wcompr (csrdtoany(B))


END SUBROUTINE reordwrap

END MODULE
