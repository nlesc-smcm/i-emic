!#begindoc

MODULE m_build

!-----------------------------------------------------------------------

!     Definition of the descriptors for the different storage types of
!     of the sparse matrices stored.

!     These storage types are:

!       Identifier  Description
!       CSC         Compact Sparse Column matrix.
!       CSR         Compact Sparse Row matrix.
!       DIAtp       (Block) Diagonal matrix; the elements in the diagonal
!                   blocks are stored in column major order.
!       FMtp        Full Matrix, stored in column major order.
!       CSRD        Compact Sparse Row, Diagonal extracted, matrix.
!       SCBM        Schur Complement from Block Matrix.
!       SCD         Schur Complement, Diagonal extracted, matrix.
!                   (CSRD with restriction)

!       MLP         Multi-Level Partioned LDU matrix, with the following
!                   level types:
!         PLDU      Partition with incomplete L, D and U decomposing
!                   factors in an MLP type matrix.
!         PFFP      Last Partition with Full matrix Factoring with
!                   Pivots.
!         PSFP      Last Partition with Sparse matrix Factoring with
!                   Pivots.

!       PRC         Preconditioner

!       Not yet implemented are:
!       BDIAR       Block Diagonal, stored by Row, matrix.
!       JDS         General sparse JDS matrix
!       SYM         Symmetric matrix with pivoting

!-----------------------------------------------------------------------
!     NOT YET IMPLEMENTED:
!-----------------------------------------------------------------------

!     bdiar	    Type identifier for a Block Diagonal, stored by Row, matrix.
!     jds	    Type identifier for a general sparse JDS matrix
!     sym	    Type identifier for a SYMmetric matrix with pivoting

!-----------------------------------------------------------------------

!#enddoc

!=======================================================================
!     Common to ALL descriptors:
!=======================================================================
!
!     	typ	Type of storage:
!	n	Dimension of matrix; number of columns and/or rows.
!
!=======================================================================

#ifdef WITH_UNION

TYPE anymatrix
      INTEGER 					:: typ
      INTEGER					:: n
END TYPE

#else

#define csrmatrix  anymatrix
#define cscmatrix  anymatrix
#define fmmatrix   anymatrix
#define diamatrix  anymatrix
#define csrdmatrix anymatrix
#define scdematrix anymatrix
#define scbmmatrix anymatrix
#define prcmatrix  anymatrix
#define partmatrix anymatrix
#define mlpmatrix  anymatrix


TYPE anymatrix
      INTEGER 					:: typ
      INTEGER					:: n, nnz

#endif

!-----------------------------------------------------------------------
!
!     CSRtp         Compact Sparse Row matrix
!     ------------------------------------------------------------------
!
!       typ	Storage Type
!       n	Dimension of matrix
!       co	Location of the nonzero matrix Values.
!                      Size of segment:  beg(N+1)-beg(1)
!       beg     Location of the begin of each Row in jco- and co- vector.
!                      Size of segment:  N+1
!       jco     Location of the Column indices in each row.
!                      Size of segment:  beg(N+1)-beg(1)
!
!-----------------------------------------------------------------------

#ifdef WITH_UNION

TYPE csrmatrix ! extends anymatrix
      INTEGER 					:: typ
      INTEGER					:: n, nnz
      INTEGER, DIMENSION(:), POINTER 		:: beg
      INTEGER, DIMENSION(:), POINTER 		:: jco
      DOUBLE PRECISION, DIMENSION(:), POINTER 	:: co
END TYPE

#else
      INTEGER, DIMENSION(:), POINTER 		:: beg
      INTEGER, DIMENSION(:), POINTER 		:: jco
      DOUBLE PRECISION, DIMENSION(:), POINTER 	:: co

#endif

!-----------------------------------------------------------------------
!
!     CSCtp         Compact Sparse Column matrix
!     ------------------------------------------------------------------
!
!       typ     Storage Type
!       n	Dimension of matrix
!       co	Location of the nonzero matrix Values.
!                      Size of segment:  beg(N+1)-beg(1)
!       beg	Location of the begin of each Column in jco- and co- vector.
!                      Size of segment:  N+1
!       jco     Location of the Row indices in each column.
!                      Size of segment:  beg(N+1)-beg(1)
!
!-----------------------------------------------------------------------

#ifdef WITH_UNION

TYPE cscmatrix ! extends anymatrix
      INTEGER 					:: typ
      INTEGER					:: n, nnz
      INTEGER, DIMENSION(:), POINTER 		:: beg
      INTEGER, DIMENSION(:), POINTER 		:: jco
      DOUBLE PRECISION, DIMENSION(:), POINTER 	:: co
END TYPE

#else

#endif

!-----------------------------------------------------------------------
!
!     DIAtp:        Sparse (block-) Diagonal, stored by Column, matrix.
!                   The elements of the (full) blocks are stored by
!                   column.
!                   {  1 <= BlkSiz , MOD (N, BlkSiz) = 0  }
!     ------------------------------------------------------------------

!       typ	Storage Type
!       n	Dimension of matrix
!                      Value:            N
!       co	Location of the matrix Values.
!                      Size of segment:  N * BlkSiz
!       blksiz	Number of rows/columns in one block in block-diagonal matrix.
!                      Value:            BlkSiz
!
!-----------------------------------------------------------------------

#ifdef WITH_UNION

TYPE diamatrix ! extends anymatrix
      INTEGER 					:: typ
      INTEGER					:: n
      INTEGER					:: blksiz
      DOUBLE PRECISION, DIMENSION(:,:), POINTER :: com
END TYPE

#else

      INTEGER					:: blksiz
      DOUBLE PRECISION, DIMENSION(:,:), POINTER :: com

#endif

!-----------------------------------------------------------------------
!
!     FMtp:         Full Matrix, stored by column.
!     ------------------------------------------------------------------
!
!       typ        Storage Type
!       n        Dimension of matrix
!                      Value:            N
!       co       Location of the matrix Values.
!                      Size of segment:  N*N
!
!-----------------------------------------------------------------------

#ifdef WITH_UNION

TYPE fmmatrix ! extends anymatrix
      INTEGER 					:: typ
      INTEGER					:: n
      DOUBLE PRECISION, DIMENSION(:,:), POINTER :: com
END TYPE

#else

#endif

!-----------------------------------------------------------------------
!
!     SCDEtp:       Schur Complement, Diagonal Extracted, matrix.
!     ------------------------------------------------------------------
!
!       typ	Storage Type
!       n	Dimension of matrix
!       dia	Location of the sparse (block-)diagonal, stored by column, matrix.
!       offd	Location of the CSR matrix containing the off-diagonal elements.
!
!-----------------------------------------------------------------------

#ifdef WITH_UNION

TYPE scdematrix ! extends anymatrix
      INTEGER 					:: typ
      INTEGER					:: n
      TYPE (diamatrix), POINTER			:: dia
      TYPE (csrmatrix), POINTER			:: offd
END TYPE

#else

      TYPE (diamatrix), POINTER			:: dia
      TYPE (csrmatrix), POINTER			:: offd

#endif

!-----------------------------------------------------------------------
!
!     CSRDtp        Compact Sparse Row matrix, Diagonal extracted
!                   The diagonal matrix may a scalar or block diagonal
!                   matrix.
!     ------------------------------------------------------------------
!
!       typ	Storage Type
!       n	Dimension of matrix
!       dia	Location of the sparse (block-)diagonal, stored by column, matrix.
!       offd	Location of the CSR matrix containing the off-diagonal elements.
!       lotr	Location of the indices of the last elements of the Lower Left
!               partition in the representation of the matrix.
!               The "location" NULLINX indicates an SCD type matrix,
!               i.e. a CSRD representation without this segment.
!
!-----------------------------------------------------------------------

#ifdef WITH_UNION

TYPE csrdmatrix ! extends scdematrix
      INTEGER 					:: typ
      INTEGER					:: n
      TYPE (diamatrix), POINTER			:: dia
      TYPE (csrmatrix), POINTER			:: offd
      INTEGER, DIMENSION(:), POINTER 		:: lotr
END TYPE

#else

      INTEGER, DIMENSION(:), POINTER 		:: lotr

#endif

!-----------------------------------------------------------------------
!
!     SCBM:         Schur Complement from Block Matrix.
!                   The matrix is partitioned into 4 blocks and the
!                   left-upper block is a (block-) diagonal matrix.
!     ------------------------------------------------------------------
!
!       typ	Storage Type
!       n       Dimension of matrix
!       g       Number of columns/rows in the left upper part, A11,
!               of a block factored matrix A.
!       nschur  = n - g
!       a11D    Location of the descriptor for the
!               inverse of the block A11, inv(A11), of a block
!               factored matrix A. A11 is a sparse (block-)diagonal,
!               stored by column, matrix.
!       a12     Location of the descriptor for the
!               right upper block (submatrix), A12, of a matrix A.
!       a21     Location of the descriptor for the
!               lower left block (submatrix), A21, of a matrix A.
!       a22     Location of the descriptor for the
!               lower right block (submatrix), A22, of a matrix A.
!
!-----------------------------------------------------------------------

#ifdef WITH_UNION

TYPE scbmmatrix ! extends anymatrix
      INTEGER 					:: typ
      INTEGER					:: n
      INTEGER					:: g
      INTEGER					:: nschur
      TYPE (diamatrix), POINTER			:: a11d
      TYPE (csrmatrix), POINTER			:: a12
      TYPE (csrmatrix), POINTER			:: a21
      TYPE (csrmatrix), POINTER			:: a22
END TYPE

#else

      INTEGER					:: g
      INTEGER					:: nschur
      TYPE (diamatrix), POINTER			:: a11d
      TYPE (csrmatrix), POINTER			:: a12
      TYPE (csrmatrix), POINTER			:: a21
      TYPE (csrmatrix), POINTER			:: a22

#endif

!=======================================================================
!     PRECONDITIONER:
!=======================================================================
!
!     PRC:          Preconditioner
!     ------------------------------------------------------------------
!
!       typ	Storage Type
!       n       Dimension of matrix
!       g       Number of rows/columns in the left upper part in the
!               partitioned Red-Black reordered matrix.
!       nschur  = n - g
!       scale   Location of the vector with Scale Factors.
!                      Size of segment:  N
!       perrb   Location of the permutation vector for the Red-Black
!                   reordering.
!                      Size of segment:  N
!       aro     Location of CSR or SCBM type matrix. The (partitioned) Red-Black
!               reordered matrix of the linear system. The matrix of
!               the linear system to be solved with a CG-type method
!               is the Schur-complement of Aro11 in Aro.
!       mlp     Location of the Multi Level Partitioned LDU matrix, OR
!               NULLINX  if not yet initialised
!                      Name:             MLPrc
!
!-----------------------------------------------------------------------

#ifdef WITH_UNION

TYPE prcmatrix ! extends anymatrix
      INTEGER 					:: typ
      INTEGER					:: n
      INTEGER					:: g
      INTEGER					:: nschur
      DOUBLE PRECISION, DIMENSION(:), POINTER	:: scale
      INTEGER, DIMENSION(:), POINTER		:: perrb
      TYPE (anymatrix), POINTER			:: aro
      TYPE (mlpmatrix), POINTER			:: mlp
END TYPE

#else

      DOUBLE PRECISION, DIMENSION(:), POINTER	:: scale
      INTEGER, DIMENSION(:), POINTER		:: perrb
      TYPE (anymatrix), POINTER			:: aro
      TYPE (mlpmatrix), POINTER			:: mlp

#endif


!=======================================================================
!     MULTI-LEVEL PARTITIONED LDU MATRIX:
!=======================================================================
!
!     MLP:          Multi-Level Partitioned LDU matrix.
!                   The partitions are stored in a doubly linked list.
!                   The descriptor of an MLP matrix is the head of
!                   this list.
!                   All partitions, except the last, are type PLDU
!                   matrices.
!                   The last partition is:
!                   . of type  PFFP, a product of Full matrices, or
!                   . of type  PSFP, a product of Sparse matrices.
!
!     ------------------------------------------------------------------
!
!       typ	Storage Type
!       n      	Dimension of matrix
!       first	Location of the descriptor of the First partition.
!       last    Location of the descriptor of the Last partition.
!       perm    Location of the Preconditioner Permutation vector.
!                      Size of segment:  N
!
!     ------------------------------------------------------------------

#ifdef WITH_UNION

TYPE mlpmatrix ! extends anymatrix
      INTEGER 					:: typ
      INTEGER					:: n
      TYPE (partmatrix), POINTER 		:: last
      TYPE (partmatrix), POINTER 		:: first
      INTEGER, DIMENSION(:), POINTER 		:: perm
END TYPE

#else

      TYPE (partmatrix), POINTER 		:: last
      TYPE (partmatrix), POINTER 		:: first
      INTEGER, DIMENSION(:), POINTER 		:: perm

#endif

!=======================================================================
!     PARTITION DESCRIPTORS:
!=======================================================================

!     Common to all partition descriptors:
!     ------------------------------------------------------------------

!     In a partition of a multiple level partitioned matrix, the 2nd
!     element of descriptor is used for the "Partition" in the
!     root matrix.

!       opOff       Partition in the root matrix, i.e. the number
!                   of rows and columns already stored in the previous
!                   partitions.

!-----------------------------------------------------------------------
!     Common to the partition descriptors PFFP and PSFP:
!     ------------------------------------------------------------------

!       opPiv       Location of the indices of the pivot elements.

!-----------------------------------------------------------------------
!     Common to the partition descriptors PLDU and PSFP:
!     ------------------------------------------------------------------

!       opDia       Location of descriptor of a
!                   (block-) Diagonal matrix (DIAtp).

!-----------------------------------------------------------------------

!     PLDU:         Partition with incomplete L, D and U decomposing
!                   factors in an MLP type matrix.
!                   All the partitions in such a matrix, except the last
!                   one, are of type PLDU.
!     ------------------------------------------------------------------

!       typ        Storage Type
!       opOff       Partition the root matrix, i.e. the number
!                   of rows and columns already stored in the previous
!                   partitions.

!       opPrev      Location of descriptor of
!                   Previous partition
!       opNext      Location of descriptor of
!                   Next partition

!       opDia       Location of descriptor of
!                   (block-) Diagonal matrix (DIAtp).
!       opLTr       Location of descriptor of
!                   Lower Triangular matrix (CSC).
!       opUTr       Location of descriptor of
!                   Upper Triangular matrix (CSR).

!-----------------------------------------------------------------------

!     PFFP:         Last partition in an MLP-type matrix:
!                   Partition with Full matrix Factored with Pivots.
!     ------------------------------------------------------------------

!       typ        Storage Type
!       opOff       Partition in the root matrix, i.e. the number
!                   of rows and columns already stored in the previous
!                   partitions.

!       opPrev      Location of descriptor of
!                   Previous partition.
!       opNext      Location of descriptor of
!                   Next partition.

!       opPiv       Location of the Pivot elements.
!                      Segment:          piv
!                      Size of segment:  NP

!       opFM        Location of descriptor of the Full
!                   Matrix containing the L and U factors.
!                      Size of segment:  NP*NP

!-----------------------------------------------------------------------


!     PSFP:         Last partition in an MLP-type matrix:
!                   Partition with Sparse matrix Factored with Pivots.
!     ------------------------------------------------------------------

!       typ        Storage Type
!       opOff       Partition in the root matrix, i.e. the number
!                   of rows and columns already stored in the previous
!                   partitions.

!       opPrev      Location of descriptor of
!                   Previous partition
!       opNext      Location of descriptor of
!                   Next partition

!       opPiv       Location of the Pivot elements.
!                      Segment:          piv
!                      Size of segment:  NP

!       opDia       Location of descriptor of the
!                   Diagonal matrix.
!       opOffD      Location of descriptor of the
!                   Off-diagonal matrix (type: CSR).
!       opLNZL      Location of the column partition information, i.e. index, in jco-
!                   and co-vectors, of Last Nonzero in row of the Lower
!                   triangular matrix
!                      Segment:          lasnz
!                      Size of segment:  NP

#ifdef WITH_UNION

TYPE partmatrix ! extends anymatrix
      INTEGER 					:: typ
      INTEGER					:: n
      TYPE (partmatrix), POINTER 		:: prev
      TYPE (partmatrix), POINTER 		:: next
      INTEGER 					:: off
      INTEGER, DIMENSION(:), POINTER 		:: piv
      TYPE (diamatrix), POINTER 		:: dia
      TYPE (fmmatrix), POINTER 			:: fm
      TYPE (cscmatrix), POINTER 		:: ltr
      TYPE (csrmatrix), POINTER 		:: offd
      INTEGER, DIMENSION(:), POINTER		:: lnzl
      TYPE (csrmatrix), POINTER 		:: utr
END TYPE

#else

      TYPE (partmatrix), POINTER 		:: prev
      TYPE (partmatrix), POINTER 		:: next
      INTEGER 					:: off
      INTEGER, DIMENSION(:), POINTER 		:: piv
      TYPE (fmmatrix), POINTER 			:: fm
      TYPE (cscmatrix), POINTER 		:: ltr
      INTEGER, DIMENSION(:), POINTER		:: lnzl
      TYPE (csrmatrix), POINTER 		:: utr
END TYPE

#endif


#ifdef WITH_UNION

TYPE dirty
  UNION
    MAP
      TYPE (anymatrix), POINTER :: any
    END MAP
    MAP
      TYPE (csrmatrix), POINTER :: csr
    END MAP
    MAP
      TYPE (cscmatrix), POINTER :: csc
    END MAP
    MAP
      TYPE (fmmatrix), POINTER :: fm
    END MAP
    MAP
      TYPE (diamatrix), POINTER :: dia
    END MAP
    MAP
      TYPE (csrdmatrix), POINTER :: csrd
    END MAP
    MAP
      TYPE (scdematrix), POINTER :: scde
    END MAP
    MAP
      TYPE (scbmmatrix), POINTER :: scbm
    END MAP
    MAP
      TYPE (prcmatrix), POINTER :: prc
    END MAP
    MAP
      TYPE (partmatrix), POINTER :: part
    END MAP
    MAP
      TYPE (mlpmatrix), POINTER :: mlp
    END MAP
  END UNION
END TYPE

#endif

INTEGER, PARAMETER :: csctp   	= -11 
INTEGER, PARAMETER :: csrtp   	= -12 
INTEGER, PARAMETER :: diatp 	= -13 
INTEGER, PARAMETER :: fmtp 	= -14 
INTEGER, PARAMETER :: csrdtp 	= -21 
INTEGER, PARAMETER :: scbmtp 	= -22 
INTEGER, PARAMETER :: scdetp 	= -23 
INTEGER, PARAMETER :: symtp 	= -26
INTEGER, PARAMETER :: jdstp 	= -27
INTEGER, PARAMETER :: bdiartp 	= -28
INTEGER, PARAMETER :: mlptp 	= -31 
INTEGER, PARAMETER :: pldutp 	= -32 
INTEGER, PARAMETER :: pffptp 	= -33 
INTEGER, PARAMETER :: psfptp 	= -34 
INTEGER, PARAMETER :: prctp 	= -41 




CONTAINS


!================================================================================================================================================

INTEGER FUNCTION cscnnz(x)

TYPE (cscmatrix), POINTER 	:: x

!     Computes the number on non zeros for a CSC matrix from x%beg

!     Arguments:
!     ==========
!     x      i   Location of the Matrix descriptor.
!     cscnnz    o   The number of nonzeros in the matrix indicated by 'x'.

  cscnnz = SUM( x%beg(2:x%n+1) - x%beg(1:x%n) )
  
END FUNCTION cscnnz

!================================================================================================================================================

INTEGER FUNCTION csrnnz(x)

TYPE (csrmatrix), POINTER 	:: x

!     Computes the number on non zeros for a CSR matrix from x%beg

!     Arguments:
!     ==========
!     x      i   Location of the Matrix descriptor.
!     csrnnz    o   The number of nonzeros in the matrix indicated by 'x'.

!     Local Parameters:
!     =================

  csrnnz = SUM( x%beg(2:x%n+1) - x%beg(1:x%n) )
  
END FUNCTION csrnnz

!================================================================================================================================================

INTEGER FUNCTION nnz(x)

USE m_dump

TYPE (anymatrix), POINTER 	:: x

!     Computes the number on non zeros for a CSC or CSR matrix from x%beg

!     Arguments:
!     ==========
!     x      i   Location of the Matrix descriptor.
!     nnz    o   The number of nonzeros in the matrix indicated by 'x'.

!     Local Parameters:
!     =================

CHARACTER (LEN=*), PARAMETER :: rounam = 'nnz'

!     Local Variables:
!     ================

INTEGER 			:: ityp
TYPE (cscmatrix), POINTER      	:: xcsc
TYPE (csrmatrix), POINTER      	:: xcsr

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)', 'Entry:', rounam
#endif

ityp = x%typ
  
IF (ityp == csctp) THEN
  xcsc => anytocsc(x)
  nnz   = cscnnz(xcsc)
ELSE IF (ityp == csrtp) THEN
  xcsr => anytocsr(x)
  nnz   = csrnnz(xcsr)
ELSE
  PRINT '(A, 2X, A, /, 3X, A, I11, 3X, A)' , 'Internal error in', rounam,  &
      'Storage type', ityp, 'not implemented!'
  CALL dump(__FILE__,__LINE__,'Not implementd')
END IF
  
END FUNCTION nnz

!================================================================================================================================================

#ifdef WITH_UNION

FUNCTION anytodia(x)

USE m_dump

TYPE (diamatrix), POINTER :: anytodia
TYPE (anymatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= diatp) CALL dump(__FILE__,__LINE__,'Illegal typecast')
y%any    => x
anytodia => y%dia

END FUNCTION

#else

FUNCTION anytodia(x)

USE m_dump

TYPE (diamatrix), POINTER :: anytodia
TYPE (anymatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= diatp) CALL dump(__FILE__,__LINE__,'Illegal typecast')
anytodia => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION anytofm(x)

USE m_dump

TYPE (fmmatrix), POINTER :: anytofm
TYPE (anymatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= fmtp) CALL dump(__FILE__,__LINE__,'Illegal typecast')
y%any   => x
anytofm => y%fm

END FUNCTION

#else

FUNCTION anytofm(x)

USE m_dump

TYPE (fmmatrix), POINTER :: anytofm
TYPE (anymatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= fmtp) CALL dump(__FILE__,__LINE__,'Illegal typecast')
anytofm => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION anytoprc(x)

USE m_dump

TYPE (prcmatrix), POINTER :: anytoprc
TYPE (anymatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= prctp) CALL dump(__FILE__,__LINE__,'Illegal typecast')
y%any    => x
anytoprc => y%prc

END FUNCTION

#else

FUNCTION anytoprc(x)

USE m_dump

TYPE (prcmatrix), POINTER :: anytoprc
TYPE (anymatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= prctp) CALL dump(__FILE__,__LINE__,'Illegal typecast')
anytoprc => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION anytocsc(x)

USE m_dump

TYPE (cscmatrix), POINTER :: anytocsc
TYPE (anymatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= csctp) CALL dump(__FILE__,__LINE__,'Illegal typecast')
y%any    => x
anytocsc => y%csc

END FUNCTION

#else

FUNCTION anytocsc(x)

USE m_dump

TYPE (cscmatrix), POINTER :: anytocsc
TYPE (anymatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= csctp) CALL dump(__FILE__,__LINE__,'Illegal typecast')
anytocsc => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION anytocsr(x)

USE m_dump

TYPE (csrmatrix), POINTER :: anytocsr
TYPE (anymatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= csrtp)  CALL dump(__FILE__,__LINE__,'Illegal typecast')
y%any    => x
anytocsr => y%csr

END FUNCTION

#else

FUNCTION anytocsr(x)

USE m_dump

TYPE (csrmatrix), POINTER :: anytocsr
TYPE (anymatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= csrtp)  CALL dump(__FILE__,__LINE__,'Illegal typecast')
anytocsr => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION anytocsrd(x)

USE m_dump

TYPE (csrdmatrix), POINTER :: anytocsrd
TYPE (anymatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= csrdtp) CALL dump(__FILE__,__LINE__,'Illegal typecast')
y%any    => x
anytocsrd => y%csrd

END FUNCTION

#else

FUNCTION anytocsrd(x)

USE m_dump

TYPE (csrdmatrix), POINTER :: anytocsrd
TYPE (anymatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= csrdtp) CALL dump(__FILE__,__LINE__,'Illegal typecast')
anytocsrd => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION anytoscbm(x)

USE m_dump

TYPE (scbmmatrix), POINTER :: anytoscbm
TYPE (anymatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= scbmtp) CALL dump(__FILE__,__LINE__,'Illegal typecast')
y%any    => x
anytoscbm => y%scbm

END FUNCTION

#else

FUNCTION anytoscbm(x)

USE m_dump

TYPE (scbmmatrix), POINTER :: anytoscbm
TYPE (anymatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= scbmtp) CALL dump(__FILE__,__LINE__,'Illegal typecast')
anytoscbm => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION anytoscde(x)

USE m_dump

TYPE (scdematrix), POINTER :: anytoscde
TYPE (anymatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer x in anytodia')
#endif

IF (x%typ /= scdetp .AND. x%typ /= csrdtp ) CALL dump(__FILE__,__LINE__,'Illegal typecast from anymatrix to scdematrix')
y%any    => x
anytoscde => y%scde

END FUNCTION

#else

FUNCTION anytoscde(x)

USE m_dump

TYPE (scdematrix), POINTER :: anytoscde
TYPE (anymatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer x in anytodia')
#endif

IF (x%typ /= scdetp .AND. x%typ /= csrdtp ) CALL dump(__FILE__,__LINE__,'Illegal typecast from anymatrix to scdematrix')
anytoscde => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION anytomlp(x)

USE m_dump

TYPE (mlpmatrix), POINTER :: anytomlp
TYPE (anymatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= mlptp) CALL dump(__FILE__,__LINE__,'Illegal typecast')
y%any    => x
anytomlp => y%mlp

END FUNCTION

#else

FUNCTION anytomlp(x)

USE m_dump

TYPE (mlpmatrix), POINTER :: anytomlp
TYPE (anymatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= mlptp) CALL dump(__FILE__,__LINE__,'Illegal typecast')
anytomlp => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION anytopart(x)

USE m_dump

TYPE (partmatrix), POINTER :: anytopart
TYPE (anymatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= pldutp  .AND. x%typ /= pffptp  .AND. x%typ /= psfptp ) CALL dump(__FILE__,__LINE__,'Illegal typecast')
y%any    => x
anytopart => y%part

END FUNCTION

#else

FUNCTION anytopart(x)

USE m_dump

TYPE (partmatrix), POINTER :: anytopart
TYPE (anymatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

IF (x%typ /= pldutp  .AND. x%typ /= pffptp  .AND. x%typ /= psfptp ) CALL dump(__FILE__,__LINE__,'Illegal typecast')

anytopart => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION fmtoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: fmtoany
TYPE (fmmatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

y%fm    => x
fmtoany => y%any

END FUNCTION

#else

FUNCTION fmtoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: fmtoany
TYPE (fmmatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

fmtoany => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION csctoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: csctoany
TYPE (cscmatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

y%csc    => x
csctoany => y%any

END FUNCTION

#else

FUNCTION csctoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: csctoany
TYPE (cscmatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

csctoany => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION csrtoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: csrtoany
TYPE (csrmatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

y%csr    => x
csrtoany => y%any

END FUNCTION

#else

FUNCTION csrtoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: csrtoany
TYPE (csrmatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

csrtoany => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION diatoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: diatoany
TYPE (diamatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

y%dia    => x
diatoany => y%any

END FUNCTION

#else

FUNCTION diatoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: diatoany
TYPE (diamatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

diatoany => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION mlptoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: mlptoany
TYPE (mlpmatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

y%mlp    => x
mlptoany => y%any

END FUNCTION

#else

FUNCTION mlptoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: mlptoany
TYPE (mlpmatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

mlptoany => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION csrdtoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: csrdtoany
TYPE (csrdmatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

y%csrd    => x
csrdtoany => y%any

END FUNCTION

#else

FUNCTION csrdtoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: csrdtoany
TYPE (csrdmatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

csrdtoany => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION prctoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: prctoany
TYPE (prcmatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

y%prc    => x
prctoany => y%any

END FUNCTION

#else

FUNCTION prctoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: prctoany
TYPE (prcmatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

prctoany => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION scdetoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: scdetoany
TYPE (scdematrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

y%scde    => x
scdetoany => y%any

END FUNCTION

#else

FUNCTION scdetoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: scdetoany
TYPE (scdematrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

scdetoany => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION scbmtoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: scbmtoany
TYPE (scbmmatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

y%scbm    => x
scbmtoany => y%any

END FUNCTION

#else

FUNCTION scbmtoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: scbmtoany
TYPE (scbmmatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

scbmtoany => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION parttoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: parttoany
TYPE (partmatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

y%part    => x
parttoany => y%any

END FUNCTION

#else

FUNCTION parttoany(x)

USE m_dump

TYPE (anymatrix), POINTER :: parttoany
TYPE (partmatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

parttoany => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION csrtocsc(x)

USE m_dump

TYPE (cscmatrix), POINTER :: csrtocsc
TYPE (csrmatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

y%csr => x
csrtocsc => y%csc

END FUNCTION

#else

FUNCTION csrtocsc(x)

USE m_dump

TYPE (cscmatrix), POINTER :: csrtocsc
TYPE (csrmatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

csrtocsc => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION csctocsr(x)

USE m_dump

TYPE (csrmatrix), POINTER :: csctocsr
TYPE (cscmatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

y%csc => x
csctocsr => y%csr

END FUNCTION

#else

FUNCTION csctocsr(x)

USE m_dump

TYPE (csrmatrix), POINTER :: csctocsr
TYPE (cscmatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

csctocsr => x

END FUNCTION

#endif



#ifdef WITH_UNION

FUNCTION csrdtoscde(x)

USE m_dump

TYPE (scdematrix), POINTER :: csrdtoscde
TYPE (csrdmatrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

y%csrd => x
csrdtoscde => y%scde

END FUNCTION

#else

FUNCTION csrdtoscde(x)

USE m_dump

TYPE (scdematrix), POINTER :: csrdtoscde
TYPE (csrdmatrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

csrdtoscde => x

END FUNCTION

#endif



#ifdef WITH_UNION

!================================================================================================================================================

FUNCTION scdetocsrd(x)

USE m_dump

TYPE (csrdmatrix), POINTER :: scdetocsrd
TYPE (scdematrix), POINTER :: x
TYPE (dirty)	    	  :: y

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

y%scde => x
scdetocsrd => y%csrd

END FUNCTION

#else

FUNCTION scdetocsrd(x)

USE m_dump

TYPE (csrdmatrix), POINTER :: scdetocsrd
TYPE (scdematrix), POINTER :: x

#ifdef DEBUG
  IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')
#endif

scdetocsrd => x

END FUNCTION

#endif




END MODULE
