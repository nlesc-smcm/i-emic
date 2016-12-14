!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif


MODULE m_wrtvsm

CONTAINS

SUBROUTINE wrtvsm (matnm, n, a, ier)

USE m_fstrlen
USE m_build
USE m_getunit

CHARACTER (LEN=*)		, INTENT(IN)		:: matnm
INTEGER				, INTENT(IN)            :: n
TYPE (csrmatrix)		, POINTER		:: a
INTEGER				, INTENT(OUT)           :: ier

!     Writes a DOUBLE PRECISION sparse square matrix A in Compressed
!     Sparse Row (CSR) format to the 3 ASCII files  <matnm>.beg ,
!     <matnm>.jco and <matnm>.co
!     Each of these files will contain one number per line!
!     The structure of the matrix can than be visualized with the
!     command  vsm <matnm>

!     Arguments:
!     ==========
!     matnm    i   Name of matrix A
!     N        i   The order of the matrix A
!     a%beg     i   a%beg(i) = index in array a%jco or a%co for the first
!                  nonzero element A(i,j) in row i  (1 <= i <= N+1)
!                  a%beg(N+1) = index of last nonzero element + 1
!     a%jco     i   a%jco(v) = column number of v-th nonzero element
!     a%co      i   a%co(v)  = value of v-th nonzero element
!                  (1 <= v <= a%beg(N+1) - 1
!     ier      o   The error code,
!                  =   0   No error
!                  =  -3   Error during closing, opening or rewinding
!                          the file <matnm>.beg
!                  =  -5   Error during writing <matnm>.beg
!                  =  -6   Error during closing, opening or rewinding
!                          the file <matnm>.jco
!                  =  -8   Error during writing <matnm>.jco
!                  =  -9   Error during closing, opening or rewinding
!                          the file <matnm>.co
!                  = -11   Error during writing <matnm>.co
!                  = -12   No logical unit available to write on an
!                          output file.

!#enddoc

!     1995-08-14  ,  Auke van der Ploeg  (original name  writeA)
!     2003-03-05  ,  Last update (Doeke de Vries)

!     Local parameter:
!     ================

INTEGER, PARAMETER :: maxfnmlen = 80

!     Local variables:
!     ================
INTEGER :: i, nmlen, unitnr, v, annz

CHARACTER (LEN=maxfnmlen) :: filenm

nmlen  = MIN (fstrlen(matnm), maxfnmlen-4)

!     Get an unopened logical unit number in UNITNR
CALL getunit (unitnr, ier)
IF ( ier /= 0 )  GO TO 988

!     Structure of the  <matnm>.beg files:
!     N+1  integers         index in jco... and co... of the first
!                           element of the row segment

filenm(1:nmlen+4) = matnm(1:nmlen) // '.beg'
OPEN (unitnr, FILE = filenm(1:nmlen+4), STATUS = 'UNKNOWN', ERR = 997)
REWIND (unitnr, ERR = 997)
WRITE (unitnr, 810, ERR = 995)  (a%beg(i), i = 1, n+1)
CLOSE (unitnr, ERR = 997)

!     Structure of the <matnm>.jco files:
!     A%NNZ  integers         the column numbers

annz = csrnnz(a)

filenm(1:nmlen+4) = matnm(1:nmlen) // '.jco'
OPEN (unitnr, FILE = filenm(1:nmlen+4), STATUS = 'UNKNOWN', ERR = 994)
REWIND (unitnr, ERR = 994)
WRITE (unitnr, 810, ERR = 992)  (a%jco(v), v = 1, annz )
CLOSE (unitnr, ERR = 994)

!     Structure of the co<matnm>.mat files:
!     A%NNZ  double precision numbers  the actual values

filenm(1:nmlen+3) = matnm(1:nmlen) // '.co'
OPEN (unitnr, FILE = matnm(1:nmlen+3), STATUS = 'UNKNOWN', ERR = 991)
REWIND (unitnr, ERR = 991)
WRITE (unitnr, 820, ERR = 989)  (a%co(v), v = 1, annz )
CLOSE (unitnr, ERR = 991)

!     Normal return:
RETURN

810  FORMAT ( i8 )
820  FORMAT ( e22.16 )

!     Error return:

988  ier = ier - 1

989  ier = ier - 2
991  ier = ier - 1

992  ier = ier - 2
994  ier = ier - 1

995  ier = ier - 2
997  ier = ier - 3

!     End of  wrtvsm
END SUBROUTINE wrtvsm

END MODULE