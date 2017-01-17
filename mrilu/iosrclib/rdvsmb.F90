!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

MODULE m_rdvsmb

CONTAINS

SUBROUTINE rdvsmb (matnm, a, ier, errinf)

USE m_fstrlen
USE m_build
USE m_getunit
USE m_wacsr

CHARACTER (LEN=*)			, INTENT(IN)		:: matnm
TYPE (csrmatrix)			, POINTER		:: a
INTEGER					, INTENT(OUT)           :: ier
INTEGER					, INTENT(OUT)           :: errinf

!     Reads the sparse matrix, A, given in old VSM format from the
!     3 ASCII files  <matnm>.beg , <matnm>.jco and <matnm>.co
!     and adds 1 to all entries of the beg array.
!     Each of these 3 files contains one number per line!
!     The nonzero elements of the sparse matrix will be stored, in
!     CSR format, in the 4 parameters
!        A%N, a%beg, a%jco and a%co

!     Arguments:
!     ==========
!     matnm    i   Common prefix of the 3 file names.
!     A%N        o   The order of the matrix A
!     a%beg     o   a%beg(i) = index in arrays 'a%jco' and 'a%co' for the
!                  first nonzero element A(i,j) in row i (1<=i<=A%N+1).
!                  a%beg(A%N+1) = index of last nonzero element + 1.
!                  The actual number of nonzero elements in A is
!                  a%beg(A%N+1) - a%beg(1).
!     a%jco     o   a%jco(v) = column number of v-th nonzero element
!     a%co      o   a%co(v)  = value of v-th nonzero element
!                            (1 <= v <= a%beg(A%N+1)-1)
!     ier      o   The error code,
!                  =  0   No error
!                  = -1   Order of matrix > UBOUND(a%beg,1)-1, the actual order is
!                         returned in 'A%N' and in 'errinf'.
!                  = -2   Number of nonzero elements > UBOUND(a%jco) or
! 			  number of nonzero elements > UBOUND(a%co), the
!                         actual number of nonzero elements is returned
!                         in 'errinf'.
!                  = -3   Error during opening, rewinding or closing
!                         the file <matnm>.beg
!                  = -4   Unexpected End of File reading <matnm>.beg
!                  = -5   Error during reading <matnm>.beg
!                  = -6   Error during opening, rewinding or closing
!                         the file <matnm>.jco
!                  = -7   Unexpected End of File reading <matnm>.jco
!                  = -8   Error during reading <matnm>.jco
!                  = -9   Error during opening, rewinding or closing
!                         the file <matnm>.co
!                  =-10   Unexpected End of File reading <matnm>.co
!                  =-11   Error during reading <matnm>.co
!                  =-12   No logical unit available to read from the
!                         files with the representation of the matrix.
!     errinf   o   Extra error information for 'ier' = -1 or -2 else
!                  undefined.

!#enddoc

!     1997-08-12  ,  Doeke de Vries
!     2003-03-05  ,  Last update (Doeke de Vries)

!     Local parameter:
!     ================

INTEGER, PARAMETER :: maxfnmlen = 80

!     Local variables:
!     ================
INTEGER :: dummy, nmlen, nr, unitnr, v, dummyn, i

CHARACTER (LEN=maxfnmlen) :: filenm

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'rdvsm'

PRINT '(A, X, A)' , 'Entry:', rounam
#endif

errinf = 0

nmlen  = MIN (fstrlen(matnm), maxfnmlen-4)

!     Get an unopened logical unit number in UNITNR
CALL getunit (unitnr, ier)
IF ( ier /= 0 )  GO TO 988


!     Structure of the  <matnm>.beg files:
!     N+1  integers         index in jco... and co... of the first
!                           element of the row segment

filenm(1:nmlen+4) = matnm(1:nmlen) // '.beg'
OPEN (unitnr, FILE = filenm(1:nmlen+4), STATUS = 'OLD', ERR = 997)
REWIND (unitnr, ERR = 997)

!     Count rows:
nr = 0
10   CONTINUE
READ (unitnr, *, END = 100, ERR = 995)  dummyn
nr = nr + 1
dummy=dummyn
GO TO 10
100  CONTINUE
!     End of file reached on <matnm>.beg
PRINT *,nr,dummy,dummyn
CALL wacsr (nr-1, dummy, a)
PRINT *,a%n
! Read .beg file
REWIND(unitnr)
DO i=1,nr  
  READ (unitnr, *, END = 995, ERR = 995)  a%beg(i)
  a%beg(i)=a%beg(i)+1
ENDDO
!! 200 STOP 'error reading beg array second time'
PRINT *,a%beg(a%n+1)
CLOSE (unitnr, ERR = 997)

IF ( nr <= 0 )  GO TO 996
  
  a%n   = nr - 1
  IF (a%n > UBOUND(a%beg,1)-1) THEN
    errinf = a%n
    GO TO 999
  END IF
  
  a%nnz = a%beg(a%n+1) - a%beg(1)
  IF ( ( a%nnz > UBOUND(a%jco,1) ) .OR.( a%nnz > UBOUND(a%co,1) ) )  THEN
    errinf = a%nnz
    GO TO 998
  END IF
  
  
!     Structure of the <matnm>.jco files:
!     A%NNZ  integers         the column numbers
  
  filenm(1:nmlen+4) = matnm(1:nmlen) // '.jco'
  OPEN (unitnr, FILE = filenm(1:nmlen+4), STATUS = 'OLD', ERR = 994)
  REWIND (unitnr, ERR = 994)
READ (unitnr, *, END = 993, ERR = 992)  (a%jco(v), v = 1, a%nnz)
CLOSE (unitnr, ERR = 994)


!     Structure of the <matnm>.co files:
!     A%NNZ  double precision numbers  the actual values

filenm(1:nmlen+3) = matnm(1:nmlen) // '.co'
OPEN (unitnr, FILE = filenm(1:nmlen+3), STATUS = 'OLD', ERR = 991)
REWIND (unitnr, ERR = 991)
READ (unitnr, *, END = 990, ERR = 989)  (a%co(v), v = 1, a%nnz)
CLOSE (unitnr, ERR = 991)

!     Normal return:
RETURN

!     Error return:

988  ier = ier - 1

989  ier = ier - 1
990  ier = ier - 1
991  ier = ier - 1

992  ier = ier - 1
993  ier = ier - 1
994  ier = ier - 1

995  ier = ier - 1
996  ier = ier - 1
997  ier = ier - 1

998  ier = ier - 1
999  ier = ier - 1

!     End of  rdvsmb
END SUBROUTINE rdvsmb

END MODULE