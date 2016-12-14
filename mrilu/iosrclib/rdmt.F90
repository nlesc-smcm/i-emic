!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

MODULE m_rdmt

CONTAINS

SUBROUTINE rdmt (filenm, binary, a)

USE m_dump
USE m_getunit
USE m_build

CHARACTER (LEN=*)			, INTENT(IN)		:: filenm
LOGICAL					, INTENT(IN)		:: binary
TYPE (csrmatrix)			, POINTER		:: a

!     Reads a sparse matrix A in Compressed Sparse Row (CSR) format
!     from the ASCII file 'filenm'.  The representation of the matrix
!     will be stored, rowwise, in the 4 parameters
!        A%N, a%beg, a%jco and a%co
!     Numbers in the file should be separated by at least one space or
!     end-of-line!
!     The order of the numbers in the file are
!        integer:                the order of the matrix N
!        End-of-line
!        N+1 * integer:          indices indicating the beginning of a row
!        End-of-line
!        NNZ * integer:          column numbers of nonzero values
!        End-of-line
!        NNZ * double precision: nonzero values

!     Arguments:
!     ==========
!     filenm   	i   Filename
!     A         o   Input: undefined
!                   Output: a CSR-matrix
!     A%N      	o   The order of the matrix A
!     a%beg     o   a%beg(i) = index in arrays 'a%jco' and 'a%co' for the
!                   first nonzero element A(i,j) in row i (1<=i<=A%N+1).
!                   a%beg(A%N+1) = index of last nonzero element + 1.
!                   The actual number of nonzero elements in A is
!                  a%nnz.
!     a%jco     o   a%jco(v) = column number of v-th nonzero element
!     a%co      o   a%co(v)  = value of v-th nonzero element
!                            (1 <= v <= a%nnz)

!#enddoc

!     1997-09-18  ,  Doeke de Vries
!     2003-03-05  ,  Last update (Doeke de Vries)

!     Local variables:
!     ================

INTEGER :: i, unitnr, v, ier, kode

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'rdmt'

PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Get an unopened logical unit number in UNITNR
CALL getunit (unitnr, ier)
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Error during reading filenm')


!     Structure of the file:
!     1    * integer                 Actual number of rows: A%N
!     A%N+1  * integer
!     A%NNZ  * integer                 Number of nonzeros: A%NNZ
!     A%NNZ  * double precision

IF (binary) THEN
  OPEN (unitnr, FILE = filenm, STATUS = 'OLD', FORM = 'UNFORMATTED', ERR = 997)
ELSE
  OPEN (unitnr, FILE = filenm, STATUS = 'OLD', ERR = 997)
END IF

REWIND (UNIT=unitnr, ERR = 997)

!     Request segment for the matrix descriptor:
ALLOCATE( a, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
a%n   = 0
a%nnz = 0
a%typ 	= csrtp

IF (binary) THEN
  READ (UNIT=unitnr, END = 996, ERR = 995) a%n
ELSE
  READ (UNIT=unitnr, FMT=*, END = 996) a%n
ENDIF

ALLOCATE( a%beg(1:a%n+1), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

IF (binary) THEN
  READ (UNIT=unitnr, END = 996, ERR = 995)  (a%beg(i), i = 1, a%n+1)
ELSE
  READ (UNIT=unitnr, FMT=*, END = 996)  (a%beg(i), i = 1, a%n+1)
ENDIF

a%nnz = a%beg(a%n+1) - a%beg(1)


ALLOCATE( a%jco(1:a%nnz), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( a%co(1:a%nnz), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

IF (binary) THEN 
  READ (UNIT=unitnr, END = 996, ERR = 995)  (a%jco(v), v = 1, a%nnz)
  READ (UNIT=unitnr, END = 996, ERR = 995)  (a%co(v), v = 1, a%nnz)
ELSE
  READ (UNIT=unitnr, FMT=*, END = 996 )  (a%jco(v), v = 1, a%nnz)
  READ (UNIT=unitnr, FMT=*, END = 996 )  (a%co(v), v = 1, a%nnz)
ENDIF

CLOSE (unitnr, ERR = 997)

!     Normal return:
RETURN

!     Error return:
995  CALL dump(__FILE__,__LINE__,'No logical unit available to read from filenm')
996  CALL dump(__FILE__,__LINE__,'Unexpected End Of File encountered')
997  CALL dump(__FILE__,__LINE__,'Error during closing, opening or rewinding file filenm')

!     End of  rdmt
END SUBROUTINE rdmt

END MODULE
