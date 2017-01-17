
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define cscmatrix  anymatrix

#endif


MODULE m_readsherman

CONTAINS

SUBROUTINE readsherman (maxneq, maxnz, n, a, rhs, ier)
 
USE m_build
USE m_wacsc
USE m_wfree
USE m_rdfilenm
USE m_fstrlen
USE m_getunit
USE m_guerrmsg
USE m_ioerrmsg
USE m_csr2csc

INTEGER					, INTENT(IN)		:: maxneq
INTEGER					, INTENT(IN)            :: maxnz
INTEGER					, INTENT(OUT)           :: n
TYPE (csrmatrix)			, POINTER		:: a
DOUBLE PRECISION, DIMENSION(maxneq)	, INTENT(OUT)         	:: rhs
INTEGER					, INTENT(OUT)           :: ier

!     Reads sparse matrix A and right-hand side of a system of linear
!     equations from ASCII files.

!     All 5 SHERMAN problem matrices with their right-hand sides
!     are stored in one large file.

!     The nonzero elements of the sparse matrix will be stored,
!  !! columnwise (!!!!) in the 4 parameters
!        [N, a%beg, a%jco, a%co]

!     Arguments:
!     ==========
!     MaxNeq  i   Maximum number of equations (order of the matrix)
!     MaxNZ   i   Maximum number of nonzero elements in matrix
!     N       o   The order of the matrix
!     beg     o   beg(i) = index in array jco or co for the first
!                 nonzero element A(i,j) in column i  (1 <= i <= N+1)
!                 beg(N+1) = index of last nonzero element + 1
!     jco     o   jco(Nz) = row number of Nz-th nonzero element
!     co      o   co(Nz)  = value of Nz-th nonzero element
!                 (1 <= Nz <= beg(N+1) - 1
!     Rhs     o   Rhs(i)   = value of i-th component of right-hand side
!                 (1 <= i <= N)


CHARACTER (LEN=*), PARAMETER :: rounam = 'readsherman'

!     Local Variables:
!     ================
INTEGER 				:: i, irow, nz, probty, iprob, unitnr, errinf
CHARACTER (LEN=12) 			:: filnam
CHARACTER (LEN=100) 			:: st
INTEGER, DIMENSION(1:5)			:: listn
TYPE (cscmatrix), POINTER          	:: xu

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

! --- Local Data Initialisation
DATA filnam /  'sherman.data' /
DATA listn  / 1000, 1080, 5005, 1104, 3312 /

ier = 0

!     Check if data can be stored into arguments:
DO i = 1, 5
  n = listn(i)
  IF (n > maxneq) THEN
    errinf = n
    GO TO 999
  END IF
END DO

!     Get a unit number for reading the file:
CALL getunit (unitnr, ier)
IF ( ier /= 0 ) THEN
  CALL guerrmsg (rounam, ier)
  GO TO 1000
END IF

OPEN  (unitnr, FILE = filnam, STATUS = 'OLD', ERR = 997)

REWIND (unitnr, ERR = 997)

!     Allocate segments to store the matrices:
CALL wacsc (n, maxnz, xu)

! === Get user input

100 CONTINUE
PRINT '(/, A)', 'Which Sherman problem (1..5)?'
READ (*, *, ERR=100) probty

! === Get to the problem number data

DO iprob = 1,probty
  
!   +++    Read the header
  READ (unitnr, '(A100)', END = 996, ERR = 995) st
  IF (iprob == probty) PRINT '(A)', st

  DO i = 1, 4
    READ (unitnr, '(A100)', END = 996, ERR = 995) st
  END DO

!          Read the data:
  n = listn(iprob)
  READ (unitnr, *, END = 996, ERR=995) (xu%beg(irow), irow = 1, n+1)

  a%nnz = xu%beg(n+1)-xu%beg(1)
  IF (a%nnz > UBOUND(xu%jco,1)) THEN
    errinf = a%nnz
    GO TO 998
  END IF

  READ (unitnr, *, END = 996, ERR = 995) (xu%jco(nz), nz = 1, a%nnz)
  READ (unitnr, *, END = 996, ERR = 995) (xu%co(nz) , nz = 1, a%nnz)

  READ (unitnr, *, END = 996, ERR = 995) (rhs(irow), irow = 1, n)

  IF (iprob == probty) PRINT *, 'N is ',n,' number of nonzero entries is ',a%nnz
END DO

CLOSE (unitnr, ERR = 997)

!     Convert matrix from CSC format to CSR format
CALL csr2csc (n, n, 1, xu, a)

CALL cscfree(xu)

!     Normal return:
RETURN


!     Fatal Error, issue error message and return:
995  ier = ier - 1
996  ier = ier - 1
997  ier = ier - 1
998  ier = ier - 1
999  ier = ier - 1

CALL ioerrmsg (filnam, rounam, 'nonsym', ier, errinf)

!     Error Return, force value ier /= 0:
1000 CONTINUE
ier = -1

!     End of  readsherman
END SUBROUTINE readsherman

END MODULE