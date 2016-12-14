!#begindoc

#ifndef WITH_UNION

#define scdematrix  anymatrix
#define prcmatrix  anymatrix
#define partmatrix  anymatrix
#define csrdmatrix  anymatrix

#endif

MODULE m_bepnum

CONTAINS
 
SUBROUTINE bepnum (BlkSiz, S, x)

USE m_build
USE m_dump
USE m_wappnd
USE m_wennz
USE m_wfree
USE m_csflpar
USE m_csslpar
USE m_glbpars
USE m_prcpars
USE m_vispars
USE m_xfnminr
USE m_fstrlen
USE m_invbdia
USE m_lumpspace
USE m_permpr
USE m_reordwrap
USE m_schurcmpl
USE m_stopldu
USE m_visscd

INTEGER			, INTENT(IN)    :: BlkSiz
TYPE (scdematrix)	, POINTER	:: S
TYPE (prcmatrix)	, POINTER	:: x

!     Computes the MLP type preconditioner (sub-)matrix, P, for the
!     (X%NSCHUR)x(X%NSCHUR) (sub-)matrix  S.
!     S is the first Schur-complement after the exact elimination of
!     'X%G' of the original 'N' unknowns (X%G > 0), or
!     the reordered matrix of the original linear system (X%G = 0).
!     On entry the matrix  S  is stored. The storage for  S
!     is released on exit.
!     The matrix  P  will be stored; the
!     location of the descriptor is 'x%mlp'.

!     Arguments:
!     ==========
!     X%N      i   Number of rows/columns in the original matrix A.
!     X%G        i   Number of rows/columns in the left upper part of A,
!                  corresponding with the eliminated unknowns.
!     BlkSiz   i   Number of rows/columns in a diagonal block.
!                  Both 'X%N' and 'X%G' should be an integer multiple of
!                  'BlkSiz'.
!     S      io  In:  Location of descriptor of the
!                       matrix S, the Schur-complement of  A_11  in the
!                       reordered matrix  A.
!                  Out: No location.
!                       (The storage for S has been released.)
!     x%mlp    io  Location of descriptor of the Multi Level
!                  Partitioned preconditioner matrix P.
!                  In:  descriptor of MLP type matrix with zero
!                       partitions.
!                  Out: descriptor of Preconditioner matrix P,
!                       stored in MLP format.
!
!#enddoc

!     Local Variables:
!     ================
!     density  Density of the last constructed Schur-complement.
!     MaxLas   Maximum number of rows/columns in Schur-complement that
!              enforces an exact decomposition of this matrix.
!     NEqDon   Number of rows/columns of the  X%NSCHUR  that have been
!              entered in the Multilevel Partitioned LDU matrix  P.
!     NPart    Number of Partitions stored in Multilevel Partitioned LDU
!              matrix [x%mlp].
!     X%nschur   Number of rows/columns in the current Schur-complement
!              matrix that has to be factored.
!              The size of the first Schur-complement is:  X%NSCHUR.
!     Nupp     Number of rows/columns in the current Schur-complement
!              which are stored into the Multilevel Partitioned LDU
!              matrix P.

!     CSpace  Location of (used) Column lump space.
!     RSpace  Location of (used) Row lump space.
!     Part   Location  of the descriptor of the actual
!              partition in the MLP preconditioner.
!     X%mlp%perm Location of renumbering permutation of the
!              constructed preconditioner.
!              ( X%mlp%perm(nieuwnummer) = nummerleesvolgorde. )

!     The parameters:
!     S      Location of the descriptor of the actual
!              Schur-complement matrix to be processed.


DOUBLE PRECISION 				:: density
INTEGER 					:: MaxLas, NEqDon, NEqNotDon, npart, nschur, nupp, ier
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: CSpace, RSpace
TYPE (partmatrix), POINTER			:: Part
CHARACTER (LEN=20) 				:: scfnm
INTEGER 					:: fnmlen
INTEGER 					:: irow
TYPE (csrdmatrix), POINTER          		:: B
INTEGER, DIMENSION(:), POINTER			:: xmlpperm

CHARACTER (LEN=*), PARAMETER 		:: rounam = 'bepnum'

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Initialisation:

npart  = 0
nupp   = x%nschur
!     Initial value Nupp ensures:  DBLE(Nupp)/DBLE(X%nschur) >= LocFrac

MaxLas = SQRT(DBLE(x%nschur))

!     Store the size of the MLP Preconditioner in descriptor:
x%mlp%n = x%nschur

!     Initialise permutation vector for the Preconditioner:
!        X%mlp%perm := 1:X%nschur
FORALL (irow = 1:x%nschur)  x%mlp%perm(irow) = irow

ALLOCATE( CSpace(1:x%nschur), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( RSpace(1:x%nschur), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

IF (clsonce) THEN
! Compute Lump Space Once.
! Set up Row- and Column Lump space for Preconditioner:
  IF (S%n /= x%Nschur) THEN
    PRINT '( /, A, I10, A, I10 /)' , 'Number rows/columns in block diagonal: ',&
    S%n, ' differs from number rows/columns in matrix: ', x%nschur
    CALL dump(__FILE__,__LINE__,'Internal error')
  END IF

  CALL lumpspace (S, CSpace, RSpace)
ELSE
! Compute Lump Space for each new Schur-complement.
! Initialise the used lump spaces with zero:
  CSpace = 0.0D0
  RSpace = 0.0D0
END IF

NEqDon = 0
NEqNotDon = x%nschur  
 
DO WHILE (.true.)

!   Density of last partition:
      
    density = DBLE(wennz (scdetoany(S))) / (DBLE(NEqNotDon)**2)

#ifdef DEBUG
    PRINT '(/, 2(A, I12  )/)', 'Done    = ', NEqDon , ' Not done = ', NEqNotDon
    PRINT '(/, 2(A, E12.6)/)', 'Density = ', density, ' Limit    = ', denslim
#endif

    IF (density   >= denslim)                  EXIT
    IF (NEqNotDon <= MaxLas)                   EXIT
    IF (DBLE(NEqNotDon)/DBLE(x%n) <= globfrac) EXIT
    IF (DBLE(nupp)/DBLE(NEqNotDon) < locfrac)  EXIT

!   Last Schur-complement is too large for an ILDU decomposition
!   and the matrix density is too low:
    
    IF (npart == 0  .ANd. visnsc > nscwrtn) THEN
!           Visualisation of Schur-complement required.
      
      nscwrtn = nscwrtn + 1
      
      CALL xfnminr ('SchurComp', 2, nscwrtn, scfnm)
      fnmlen = fstrlen(scfnm)
      
      CALL visscd ('Visualize Schur-complement with:', scfnm(1:fnmlen), S)
    END IF
    
!   Calculate a new reordering of Schur-complement, S, and
!   place the reordered matrix in  B.  This matrix, B, will be
!   stored. The global permutation, in x%mlp%perm(NEqDon+1) is maintained
!   and the (used) lump spaces are updated:

    xmlpperm => x%mlp%perm(NEqDon+1:x%nschur)
    CALL reordwrap (NEqNotDon, nupp, BlkSiz, S, B, xmlpperm, CSpace(NEqDon+1:x%nschur), RSpace(NEqDon+1:x%nschur))

!   Release workspace occupied by matrix S:

    CALL scdefree(S)
      
    IF (nupp == 0 .OR. nupp == NEqNotDon) THEN

!     The upper partition is all of the matrix;  set last
!     Schur-complement and exit the partitioning loop:
      S => csrdtoscde(B)
        
      EXIT
    END IF
      
!   Invert left upper partition of block factored matrix,  B:
      
!   Invert left upper partition of block diagonal matrix  Bd: Bd(1:Nupp,1:Nupp) := inv(Bd(1:Nupp,1:Nupp))

    CALL invbdia (nupp, B%dia)
      
!   Store partition L, D and U factors of preconditioner into

    CALL stopldu (NEqDon, NEqNotDon, nupp, B, x%mlp%perm, Part)

!   new partition:
      
!   Insert the partition 'Part' at the end of the doubly linked list in 'x%mlp':

    CALL wappnd (Part, x%mlp)
        
    IF (outlev >= 4) PRINT '(A, I3, A, 2(I7, 1X, A))', 'Step', npart + 1, ':', &
     nupp, 'of the', NEqNotDon, 'unknowns eliminated approximately.'
        
!   Calculate an approximation of the Schur-complement of the
!   submatrix  B_11  in the partitioned matrix B and
!   store this Schur-complement into S:

    CALL schurcmpl (nupp, csrdtoany(B), schtol, S)
          
    IF (visnsc > nscwrtn) THEN

!     Visualisation of Schur-complement required.
             
      nscwrtn = nscwrtn + 1
            
      CALL xfnminr ('SchurComp', 2, nscwrtn, scfnm)
      fnmlen = fstrlen(scfnm)
            
        CALL visscd ('Visualize Schur-complement with:', scfnm(1:fnmlen), S)
    END IF


!   Free workspace of matrix  B:

    CALL csrdfree(B)
          
    NEqNotDon = NEqNotDon - nupp
    NEqDon    = NEqDon    + nupp
          
    npart = npart + 1
        
!   End WHILE-loop
  END DO
      
  DEALLOCATE( CSpace, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
  DEALLOCATE( RSpace, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
      
!   Permute column numbers in the right upper blocks, U, and the row
!   numbers in the left lower blocks, L, of all PLDU type partitions,
!   in the MLP matrix  P.
!   The inverse permutation of 'X%mlp%perm' is applied.
    
!   Permute row/column numbers in off-diagonal partitions:
    CALL permpr (x%mlp)
    
    IF (vislsc) THEN
!     Visualisation of Last Schur-complement required.
      
      CALL visscd ('Visualize last Schur-complement with:', 'LastSchurComp', S)
    END IF

    IF (NEqNotDon > 0) THEN

!     Compute Last Partition and store it
!     at 'Part', and release workspace occupied by matrix S:
!      IF ( density <= sparslim .AND. NEqNotDon >= 0 ) THEN
      IF ( NEqNotDon >= 0 ) THEN

!       Compute and Store the Sparse Last Partition:
      
        CALL csslpar (NEqDon, NEqNotDon, BlkSiz, S, Part)

      ELSE

!       Compute and Store the Full Last Partition:

        CALL csflpar (NEqDon, NEqNotDon, S, Part)

      END IF

!     Append the  partition 'Part' to the end of the doubly linked
!     list in 'ixP':
    
      CALL wappnd (Part, x%mlp)
    END IF

END SUBROUTINE bepnum

END MODULE
