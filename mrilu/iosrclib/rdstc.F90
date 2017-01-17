!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

MODULE m_rdstc

CONTAINS

SUBROUTINE rdstc (filenm, binary, a, b, ier, errinf)

USE m_getunit
USE m_build

CHARACTER (LEN=*)			, INTENT(IN)		:: filenm
LOGICAL					, INTENT(IN)		:: binary
TYPE (csrmatrix)			, POINTER	        :: a
DOUBLE PRECISION, DIMENSION(1:)		, INTENT(OUT)           :: b
INTEGER					, INTENT(OUT)           :: ier
INTEGER					, INTENT(OUT)           :: errinf

!     Reads the sparse coefficient matrix, A, and right hand side, b, of
!     a system of linear equations from the ASCII file 'filenm'.
!     The system is stored in "stencil"-format.
!     The nonzero elements of the sparse matrix will be stored, rowwise,
!     in the 4 parameters
!        A%N, a%beg, a%jco and a%co
!     The right hand side will be stored in  b.

!     Arguments:
!     ==========
!     filenm   i   Filename of the input file containg the equation in
!                  the format described in the implementation part.
!     A%N        o   The actual order of the matrix A
!     a%beg     o   a%beg(i) = index in arrays 'a%jco' and 'a%co' for the
!                  first nonzero element A(i,j) in row i (1<=i<=A%N+1).
!                  a%beg(A%N+1) = index of last nonzero element + 1.
!                  The actual number of nonzero elements in A is
!                  a%nnz.
!     a%jco     o   a%jco(v) = column number of v-th nonzero element
!     a%co      o   a%co(v)  = value of v-th nonzero element
!                            (1 <= v <= a%nnz)
!     b        o   b(i)    = value of i-th component of right hand side
!                            (1 <= i <= A%N)
!     ier      o   The error code,
!                  =  0   No error.
!                  = -1   Order of matrix > UBOUND(a%beg,1)-1 or
!                         order of matrix > UBOUND(b,1), the actual order is
!                         returned in 'N' and in 'errinf'.
!                  = -2   Number of nonzero elements > UBOUND(a%jco,1) or
!                         number of nonzero elements > UBOUND(a%co,1), the
!                         estimated number of nonzero elements is
!                         returned in 'errinf'.
!                  = -3   Error during closing, opening or rewinding
!                         file 'filenm'.
!                  = -4   Unexpected End Of File encountered
!                  = -5   Error during reading 'filenm'.
!                  =-12   No logical unit available to read from
!                         'filenm'.
!                  =-13   Inconsistent file contents
!     errinf   o   Extra error information for 'ier' = -1 or -2 else
!                  undefined.

!#enddoc

!     1997-05-12  ,  Doeke de Vries
!     2003-03-05  ,  Last update (Doeke de Vries)

!     Local variables:
!     ================
INTEGER :: i, j, ni, nj, ir, jr, rowa, unitnr, annz
DOUBLE PRECISION :: cn, ce, cs, cw, cc

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'rdstc'

PRINT '(A, X, A)' , 'Entry:', rounam
#endif

errinf = 0

!     Get an unopened logical unit number in UNITNR
CALL getunit (unitnr, ier)
IF ( ier /= 0 )  GO TO 988

a%n   = 0

!     Structure of the input file:
!     integer  integer             Ni  Nj   Nr of gridpoints:  Ni x Nj
!                                  Blankline
!     Nj-2  groups, each group containing
!        Ni-2  lines, each line containing
!           2 integers   6 reals   i, j  and  Cn, Ce, Cs, Cw, Cc  and  rhs
!                                  of the stencil in grid-point (i,j).
!                                  2 <= i <= Ni-1 .AND. 2 <= j <= Nj-1
!     Cc .ne. 0  in every grid-point (i,j)

IF (binary) THEN
  OPEN (unitnr, FILE = filenm, STATUS = 'OLD', FORM = 'UNFORMATTED', ERR = 997)
ELSE
  OPEN (unitnr, FILE = filenm, STATUS = 'OLD', ERR = 997)
END IF

REWIND (unitnr, ERR = 997)
READ (unitnr, *, END = 996, ERR = 995 ) ni, nj

a%n = (ni-2)*(nj-2)
IF ( (a%n > UBOUND(a%beg,1)-1) .OR. (a%n > UBOUND(b,1)) ) THEN
  errinf = a%n
  GO TO 999
END IF

rowa         = 0
annz         = 0
a%beg(rowa+1) = annz + 1

DO   j = 2, nj-1
  DO   i = 2, ni-1
    IF ( (annz+5 > UBOUND(a%jco,1) ) .OR. (annz+5 > UBOUND(a%co,1) ) ) THEN
!            Estimation of the number of nonzeros in A:
      errinf = annz + 5 * (a%n - rowa)
      GO TO 998
    END IF

    rowa = rowa + 1
  READ (unitnr, *, END = 996, ERR = 995) ir, jr, cn, ce, cs, cw, cc, b(rowa)
  
  IF ( ir /= i  .Or. jr /= j ) GO TO 987
    
    IF ( cs /= 0 ) THEN
      annz = annz + 1
      a%jco(annz) = rowa - (ni-2)
      a%co(annz)  = cs
    END IF
    
    IF ( cw /= 0 ) THEN
      annz = annz + 1
      a%jco(annz) = rowa - 1
      a%co(annz)  = cw
    END IF
    
!         Cc .ne. 0D0   Diagonal element
    
    IF ( cc /= 0) THEN
      annz = annz + 1
      a%jco(annz) = rowa
      a%co(annz)  = cc
    ELSE
      GO TO 987
    END IF
    
    IF ( ce /= 0 ) THEN
      annz = annz + 1
      a%jco(annz) = rowa + 1
      a%co(annz)  = ce
    END IF
    
    IF ( cn /= 0 ) THEN
      annz = annz + 1
      a%jco(annz) = rowa + (ni-2)
      a%co(annz)  = cn
    END IF
    
    a%beg(rowa+1) = annz + 1
  END DO
  
END DO
a%nnz=annz

CLOSE (unitnr, ERR = 997)

!     Normal return:
RETURN

!     Error return:
987  ier = ier - 1
988  ier = ier - 7
995  ier = ier - 1
996  ier = ier - 1
997  ier = ier - 1
998  ier = ier - 1
999  ier = ier - 1

!     end of  rdstc
END SUBROUTINE rdstc

END MODULE