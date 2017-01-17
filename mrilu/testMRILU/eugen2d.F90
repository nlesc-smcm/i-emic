!      Algemenere discretisatie op rooster met stretching
!      E.F.F.Botta@math.rug.nl


!=======================================================================

#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

MODULE m_eugen2d

CONTAINS

SUBROUTINE eugen2d (probty, mint, strtch, cvx, cvy, n, x, Xex)

USE m_dump
USE m_build
USE m_wacsr
USE m_filgrd

INTEGER				, INTENT(IN)		:: probty
INTEGER				, INTENT(IN)         	:: mint
DOUBLE PRECISION		, INTENT(IN)        	:: strtch
DOUBLE PRECISION		, INTENT(IN)         	:: cvx
DOUBLE PRECISION		, INTENT(IN)         	:: cvy
INTEGER				, INTENT(OUT)           :: n
TYPE (csrmatrix)		, POINTER		:: x
DOUBLE PRECISION, DIMENSION(:)	, POINTER 		:: Xex

!     Compute the matrix and the exact solution for the convection-
!     diffusion equation:
!        - Delta U(x,y) + CVX U_x(x,y) + CVY U_y(x,y) = F(x,y)
!     on the square [0,1]x[0,1].
!     The same grid is used in the x- and y-direction.
!     The matrix will be stored in CSR format at [x].
!     The exact solution will be stored into at [Xex].

!     Arguments:
!     ==========
!     ProbTy   i   Problem type (ProbTy = 1,2,3,4,5,6)
!                  Type of boundary conditions on x=0, y=0, x=1, y=1 are
!                  indicated by D (=Dirichlet), or N (=Neumann):
!                  ProbTy=1: one_sided stretching, DDDD
!                  ProbTy=2: one_sided stretching, DDNN
!                  ProbTy=3: one_sided stretching, DDNN
!                            (Border between grid-points)
!                  ProbTy=4: two_sided stretching, DDDD
!                  ProbTy=5: two_sided stretching, NNNN
!                  ProbTy=6: two_sided stretching, NNNN
!                            (Border between grid-points)
!     Mint     i   Number of intervals in x- and y-direction.
!     Strtch   i   Grid stretching, with:  1.0D-4 <= Strtch <= 1.0D+4
!     CVX      i   Convection factor, the coefficient of  U_x.
!     CVY      i   Convection factor, the coefficient of  U_y.
!     N        o   Number of equations in linear system.
!     x      o   Location of the descriptor of the matrix stored in
!                  CSR format.
!     Xex    o   Location of the segment containing the exact
!                  solution in the grid-points.
!                  Xex(1:N)

!     Local Parameters:
!     =================
!     Maximum number of points in the i- and j-direction of the domain.

!INTEGER, PARAMETER :: nm = 1024

!     alpha = 1; %discretisatie van y' zonder diagonaalbijdrage
!     alpha = 0; %discretisatie van y' met diagonaalbijdrage

DOUBLE PRECISION, PARAMETER :: alpha = 1.0D0

CHARACTER (LEN=*), PARAMETER :: rounam = 'eugen2d'


!     Local Variables:
!     ================

INTEGER 				:: m, ier, xnnz, nm
INTEGER 				:: i, j, eqnnum, minind, maxind
LOGICAL 				:: twosid, shiftl, shiftr
DOUBLE PRECISION 			:: a, b, s
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: h, xc, cs, cw, ce, cn, s1

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif
nm=mint+2
ALLOCATE( h(1:nm), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( xc(1:nm), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( cs(1:nm), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( cw(1:nm), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( ce(1:nm), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( cn(1:nm), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( s1(1:nm), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

m = mint

IF (mint < 3) THEN
  PRINT '(A, X, A, A, /, 3X, A)', 'Fatal error in', rounam, '!',  &
      'Choose number of intervals >= 3.'
  CALL dump(__FILE__,__LINE__,'Fatal error')
END IF
IF (mint+1 >= nm) THEN
  PRINT '(A, X, A, A, /, 3X, A, I11)', 'Fatal error in', rounam, '!',  &
      'Choose number of intervals <', nm-1
  CALL dump(__FILE__,__LINE__,'Fatal error')
END IF


!     Compute the positions of the grid points, in 'xc', and the sizes
!     of the intervals in 'h'.
!     The array xc is required for generating the exact solution and
!     the right hand side of the linear system.

twosid = (probty > 3)
shiftl = (probty == 6)
shiftr = (probty == 3 .OR. probty == 6)

CALL filgrd (mint, strtch, twosid, shiftl, shiftr, xc, h)

!     De matrix wordt gevuld zodanig dat een discretisatie
!     ontstaat van een convectie-diffusie vgl. op een vierkant
!     met (deels) Dirichlet of Neumann rvw'n.
!     Er wordt gediscretiseerd vanuit de behoudswet per volume
!     zodat het diffusieve stuk een symmetrische matrix levert.
!     Vanwege de diverse randvw'n verschilt het aantal
!     x,y-intervallen m van het aantal x,y-onbekenden M !!!!

minind = 2
maxind = mint
IF (probty == 2 .OR. probty == 5) maxind = mint+1
IF (probty == 5) minind = 1

m = maxind-minind+1
n = m*m


!     Compute and store the matrix elements:

DO i = minind, maxind
  b = h(i)*h(i+1) * ( h(i+1)+h(i) )
  a = alpha*h(i)*h(i+1)
  cw(i) = (cvx*(-a-(1-alpha)*h(i+1)*h(i+1)) - 2*h(i+1))/b
  cs(i) = (cvy*(-a-(1-alpha)*h(i+1)*h(i+1)) - 2*h(i+1))/b
  ce(i) = (cvx*(+a+(1-alpha)*h(i)  *h(i)  ) - 2*h(i)  )/b
  cn(i) = (cvy*(+a+(1-alpha)*h(i)  *h(i)  ) - 2*h(i)  )/b
  
!        S1 halved to get the same matrix as from Matlab program
  s1(i) = (h(i)+h(i+1))/2
END DO

IF (probty == 2 .OR. probty == 5) THEN
  s1(maxind) = s1(maxind)/2.0D0
  cw(maxind) = cw(maxind) + ce(maxind)
  cs(maxind) = cs(maxind) + cn(maxind)
  ce(maxind) = 0.0D0
  cn(maxind) = 0.0D0
END IF

IF (probty == 5) THEN
  s1(minind) = s1(minind)/2.0D0
  ce(minind) = cw(minind) + ce(minind)
  cn(minind) = cn(minind) + cs(minind)
  cw(minind) = 0.0D0
  cs(minind) = 0.0D0
END IF

IF (probty == 3 .OR. probty == 6) THEN
  ce(maxind) = 0.0D0
  cn(maxind) = 0.0D0
END IF
IF (probty == 6) THEN
  cw(minind) = 0.0D0
  cs(minind) = 0.0D0
END IF


!     Allocate the storage for the CSR representation of the matrix:
CALL wacsr (n, m*m + 4*m*(m-1), x)

!     Request storage for the exact solution:
ALLOCATE ( Xex(1:n), STAT=ier)
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Construct the matrix:
eqnnum = 0
xnnz  = 0
x%beg(1) = xnnz + 1

DO j = minind, maxind
  DO i = minind, maxind
    eqnnum = eqnnum + 1
    
    s = s1(i)*s1(j)
    IF (j > minind) THEN
      xnnz = xnnz+1
      x%co(xnnz)  = cs(j)*s
      x%jco(xnnz) = eqnnum-m
    END IF
    IF (i > minind) THEN
      xnnz = xnnz+1
      x%co(xnnz)  = cw(i)*s
      x%jco(xnnz) = eqnnum-1
    END IF
    xnnz = xnnz+1
    x%jco(xnnz) = eqnnum
    x%co(xnnz)  = -((cw(i)+ce(i)+cs(j)+cn(j)) * s)
    IF (i < maxind) THEN
      xnnz = xnnz+1
      x%co(xnnz)  = ce(i)*s
      x%jco(xnnz) = eqnnum+1
    END IF
    IF (j < maxind) THEN
      xnnz = xnnz+1
      x%co(xnnz)  = cn(j)*s
      x%jco(xnnz) = eqnnum+m
    END IF
    
    x%beg(eqnnum+1) = xnnz+1
    
!           Generate the starting vector for the iteration:
    
    Xex(eqnnum) = xactsol (xc(i), xc(j))
  END DO
  
END DO
x%nnz = xnnz

DEALLOCATE( h, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( xc, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( cs, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( cw, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( ce, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( cn, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( s1, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

RETURN
END SUBROUTINE eugen2d

!=======================================================================

DOUBLE PRECISION FUNCTION xactsol (x, y)



DOUBLE PRECISION, INTENT(IN)         :: x
DOUBLE PRECISION, INTENT(IN)         :: y


!     Returns the value of the exact solution of the
!     convection-diffusion equation in the point (x,y).
!     To be used in the subroutine 'eugen2d'.

xactsol = (x - 0.5D0)**2 * (y - 0.5D0)**2
!      XactSol = (x*y*(1.0D0-x)*(1.0D0-y))**2*exp(x*x*y)

!     End of  XactSol
END FUNCTION xactsol

END MODULE
