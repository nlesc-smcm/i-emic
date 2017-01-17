!#begindoc

#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

MODULE m_wrdstc

CONTAINS
 
SUBROUTINE wrdstc (filenm, binary, x, b, ier)

USE m_build
USE m_wacsr
USE m_wcompr
USE m_getunit
USE m_ioerrmsg

CHARACTER (LEN=*)		, INTENT(IN)			:: filenm
LOGICAL				, INTENT(IN)			:: binary
TYPE (csrmatrix)		, POINTER       		:: x
DOUBLE PRECISION, DIMENSION(:)	, POINTER			:: b
INTEGER				, INTENT(OUT)   		:: ier

!     Reads the sparse coefficient matrix, A, and right hand side, b, of
!     a system of linear equations from the ASCII file 'filenm'. This
!     linear system is stored in "stencil"-format on the file.
!     The matrix will be stored, in CSR format, into newly allocated
!     segments. A new descriptor is filled so
!     that this matrix can be referenced through the descriptor at
!     'x'.
!     The right hand side vector will be stored in a newly allocated
!     segment at 'b'.

!     Structure of the input file:
!     integer  integer             Ni  Nj   Nr of gridpoints:  Ni x Nj
!     End-of-line
!     Blankline
!     End-of-line
!     Nj-2  groups, each group containing:
!        Ni-2  lines, each line containing:
!           2 integers   6 reals   i, j  and  Cn, Ce, Cs, Cw, Cc  and  rhs
!                                  of the stencil in grid-point (i,j).
!                                  2 <= i <= Ni-1 .AND. 2 <= j <= Nj-1
!                                  Cc .NE. 0  in every grid-point (i,j)

!     Arguments:
!     ==========
!     filenm	i   Filename of the input file containg the equation in
!                   the format described in the implementation part.
!     x   	o   Location of the matrix descriptor.
!     b      	o   Location of the right hand side vector.
!     ier      	o   Error code,
!                   =  0   No error occurred, subroutine allocated storage
!                          as described above.
!                   <  0   Some error occurred, an error message has been
!                          written to standard output.

!#enddoc

!     1997-05-12  ,  Doeke de Vries
!     2003-03-05  ,  Last update (Doeke de Vries)

!     Local variables:
!     ================

INTEGER 			:: errinf, n
INTEGER 			:: nnza, i, j, ni, nj, ir, jr, rowa, unitnr
DOUBLE PRECISION 		:: cn, ce, cs, cw, cc

CHARACTER (LEN=*), PARAMETER :: rounam = 'wrdstc'

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

errinf = 0

!     Get an unopened logical unit number in UNITNR
CALL getunit (unitnr, ier)
IF ( ier /= 0 )  GO TO 988

n   = 0

!     Structure of the input file:
!     integer  integer             Ni  Nj   Nr of gridpoints:  Ni x Nj
!     End-of-line
!     Blankline
!     End-of-line
!     Nj-2  groups, each group containing:
!        Ni-2  lines, each line containing:
!           2 integers   6 reals   i, j  and  Cn, Ce, Cs, Cw, Cc  and  rhs
!                                  of the stencil in grid-point (i,j).
!                                  2 <= i <= Ni-1 .AND. 2 <= j <= Nj-1
!                                  Cc .ne. 0  in every grid-point (i,j)

IF (binary) THEN
  OPEN (unitnr, FILE = filenm, STATUS = 'OLD', FORM = 'UNFORMATTED', ERR = 997)
ELSE
  OPEN (unitnr, FILE = filenm, STATUS = 'OLD', ERR = 997)
END IF

REWIND (unitnr, ERR = 997)

READ (unitnr, *, END = 996, ERR = 995 ) ni, nj

IF (ni <= 2  .Or. nj <= 2) THEN
  PRINT '(/, A, 2X, A, /, 3X, A, /, 3X, A, I11, 3X, A, I11)',  &
      'Fatal error occurred in', rounam, 'Illegal number of grid points:',  &
      'Ni =', ni, 'Nj =', nj
  GO TO 1000
END IF

!     The number of unknowns = the order of the matrix:
n = (ni-2)*(nj-2)

!     Allocate a segment for the right hand side vector, b:
ALLOCATE( b(1:n), STAT=ier)
IF (ier < 0) GO TO 1000

!     Allocate storage for a CSR-type matrix of order N and
!     5*N nonzeros:
CALL wacsr (n, 5*n, x)

!     Read the matrix- and the right hand side elements. Store the
!     non-zero matrix elements and the right hand side elements.

rowa = 0
nnza = 0
x%beg(rowa+1) = nnza + 1

DO   j = 2, nj-1
  DO   i = 2, ni-1
    rowa = rowa + 1
    READ (unitnr, *, END = 996, ERR = 995) ir, jr, cn, ce, cs, cw, cc, b(rowa)
  
  IF ( ir /= i  .Or. jr /= j ) GO TO 987
    
    IF ( cs /= 0 ) THEN
      nnza = nnza + 1
      x%jco(nnza) = rowa - (ni-2)
      x%co(nnza)  = cs
    END IF
    
    IF ( cw /= 0 ) THEN
      nnza = nnza + 1
      x%jco(nnza) = rowa - 1
      x%co(nnza)  = cw
    END IF
    
!         Cc .ne. 0D0   Diagonal element
    
    IF ( cc /= 0) THEN
      nnza = nnza + 1
      x%jco(nnza) = rowa
      x%co(nnza)  = cc
    ELSE
      GO TO 987
    END IF
    
    IF ( ce /= 0 ) THEN
      nnza = nnza + 1
      x%jco(nnza) = rowa + 1
      x%co(nnza)  = ce
    END IF
    
    IF ( cn /= 0 ) THEN
      nnza = nnza + 1
      x%jco(nnza) = rowa + (ni-2)
      x%co(nnza)  = cn
    END IF
    
    x%beg(rowa+1) = nnza + 1
  END DO
  
END DO

CLOSE (unitnr, ERR = 997)

!     Compress the storage of the CSR representation of the matrix
!     and adjust the required size:
CALL wcompr (csrtoany(x))

!     Normal return:
RETURN

!     Error return:
987  ier = ier - 1
988  ier = ier - 7
995  ier = ier - 1
996  ier = ier - 1
997  ier = ier - 3

errinf = 0

CALL ioerrmsg (filenm, rounam, 'UNKNOWN?', ier, errinf)


!     Error Return, force value ier < 0:
1000 CONTINUE
ier = -1

!     end of  wrdstc
END SUBROUTINE wrdstc

END MODULE