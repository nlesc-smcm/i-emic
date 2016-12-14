
#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

MODULE m_arthura

CONTAINS

SUBROUTINE arthura (ProbTy, Mint, Strtch, CVX, Source, n, x, RL)
 
USE m_build
USE m_filgrd
USE m_wacsr

INTEGER				, INTENT(IN)		:: ProbTy
INTEGER				, INTENT(IN)            :: Mint
DOUBLE PRECISION		, INTENT(IN)         	:: Strtch
DOUBLE PRECISION		, INTENT(IN)         	:: CVX
DOUBLE PRECISION		, INTENT(IN)            :: Source
INTEGER				, INTENT(OUT)           :: n
TYPE (csrmatrix)		, POINTER               :: x
DOUBLE PRECISION, DIMENSION(:)	, POINTER  		:: RL

!     This routine produces the same matrix as the matlab function
!     "discr2.m"

!     The same grid is used in the x- and y-direction.
!     The matrix will be stored in CSR format at [x].
!     The right hand side will be stored at [RL].

!     Arguments:
!     ==========
!     ProbTy   	i   Problem type (ProbTy = 1,2,3,4,5,6)
!                   Type of boundary conditions on x=0, y=0, x=1, y=1 are
!                   indicated by D (=Dirichlet), or N (=Neumann):
!                   ProbTy=1: one_sided stretching, DDDD
!                   ProbTy=2: one_sided stretching, DDNN
!                   ProbTy=3: one_sided stretching, DDNN
!                            (Border between grid-points)
!                   ProbTy=4: two_sided stretching, DDDD
!                   ProbTy=5: two_sided stretching, NNNN
!                   ProbTy=6: two_sided stretching, NNNN
!                            (Border between grid-points)
!     Mint     	i   Number of intervals in x- and y-direction.
!     Strtch   	i   Grid stretching, with:  1.0D-4 <= Strtch <= 1.0D+4
!     CVX      	i   Convection factor, the coefficient of  U_x.
!     Source   	i
!     N        	o   Number of equations in linear system.
!     x      	o   Location of the descriptor of the matrix stored in CSR format.
!     RL     	o   Location of the segment containing the right hand side of the linear system.
!                   RL(1:N)

!     Local Parameters:
!     =================
!     Maximum number of points in the i- and j-direction of the domain.

INTEGER, PARAMETER :: nm = 1024


INTEGER, PARAMETER :: minint = 2

!     alpha = 1; %discretisatie van y' zonder diagonaalbijdrage
!     alpha = 0; %discretisatie van y' met diagonaalbijdrage

DOUBLE PRECISION, PARAMETER :: alpha = 1.0D0

CHARACTER (LEN=*), PARAMETER :: rounam = 'arthura'


!     Local Variables:
!     ================

INTEGER 					:: m, ier, xnnz
INTEGER 					:: i, j, eqnnum, minind, maxind
INTEGER 					:: minprob
LOGICAL 					:: twosid, shiftl, shiftr
DOUBLE PRECISION 				:: a, b, s
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: h, xc, cw, ce, s1

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

ALLOCATE( h(1:nm), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( xc(1:nm), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( cw(1:nm), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( ce(1:nm), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( s1(1:nm), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

m = Mint

IF ((ProbTy == 1).OR.(ProbTy == 4)) THEN
  minprob = minint+1
ELSE
  minprob = minint
END IF

IF (Mint+1 >= nm) THEN
  PRINT '(A, X, A, A, /, 3X, A, I11)', 'Fatal error in', rounam, '!',  &
      'Choose number of intervals <', nm-1
  CALL dump(__FILE__,__LINE__,'Fatal erro')
END IF
IF (Mint < minprob) THEN
  PRINT '(A, X, A, A, /, 3X, A, I11)', 'Fatal error in', rounam, '!',  &
      'Choose number of intervals >=', minprob
  CALL dump(__FILE__,__LINE__,'Fatal erro')
END IF


!     Compute the positions of the grid points, in 'xc', and the sizes
!     of the intervals in 'h'.
!     The array xc is required for generating the exact solution and
!     the right hand side of the linear system.

twosid = (ProbTy > 3)
shiftl = (ProbTy == 6)
shiftr = (ProbTy == 3 .OR. ProbTy == 6)

CALL filgrd (Mint, Strtch, twosid, shiftl, shiftr, xc, h)


!     De matrix wordt gevuld zodanig dat een discretisatie
!     ontstaat van een convectie-diffusie vgl. op een vierkant
!     met (deels) Dirichlet of Neumann rvw'n.
!     Er wordt gediscretiseerd vanuit de behoudswet per volume
!     zodat het diffusieve stuk een symmetrische matrix levert.
!     Vanwege de diverse randvw'n verschilt het aantal
!     x,y-intervallen m van het aantal x,y-onbekenden M !!!!

minind = 2
maxind = Mint
IF (ProbTy == 2 .OR. ProbTy == 5) maxind = Mint+1
IF (ProbTy == 5) minind = 1

m = maxind-minind+1
n = m*m

!     Compute and store the matrix elements:

DO i = minind, maxind
  b = h(i)*h(i+1) * ( h(i+1)+h(i) )
  a = alpha*h(i)*h(i+1)
  cw(i) = (CVX*xc(i)*(-a-(1-alpha)*h(i+1)*h(i+1))-2*h(i+1))/b
  ce(i) = (CVX*xc(i)*(+a+(1-alpha)*h(i)  *h(i)  )-2*h(i)  )/b
  
!        S1 halved to get the same matrix as from Matlab program
  s1(i) = (h(i)+h(i+1))/2
END DO

IF (ProbTy == 2 .OR. ProbTy == 5) THEN
  s1(maxind) = s1(maxind)/2.0D0
  cw(maxind) = cw(maxind) + ce(maxind)
  ce(maxind) = 0.0D0
END IF

IF (ProbTy == 5) THEN
  s1(minind) = s1(minind)/2.0D0
  ce(minind) = cw(minind) + ce(minind)
  cw(minind) = 0.0D0
END IF

IF (ProbTy == 3 .OR. ProbTy == 6) THEN
  ce(maxind) = 0.0D0
END IF
IF (ProbTy == 6) THEN
  cw(minind) = 0.0D0
END IF


!     Allocate the storage for the CSR representation of the matrix:
CALL wacsr (n, m*m + 4*m*(m-1), x)

!     Request storage for the right hand side:
ALLOCATE ( RL(1:n), STAT=ier)
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
      x%co(xnnz)  = cw(j)*s
      x%jco(xnnz) = eqnnum-m
    END IF
    IF (i > minind) THEN
      xnnz = xnnz+1
      x%co(xnnz)  = cw(i)*s
      x%jco(xnnz) = eqnnum-1
    END IF
    xnnz = xnnz+1
    x%jco(xnnz) = eqnnum
    x%co(xnnz)  = -((-Source+cw(i)+ce(i)+cw(j)+ce(j)) * s)
    IF (i < maxind) THEN
      xnnz = xnnz+1
      x%co(xnnz)  = ce(i)*s
      x%jco(xnnz) = eqnnum+1
    END IF
    IF (j < maxind) THEN
      xnnz = xnnz+1
      x%co(xnnz)  = ce(j)*s
      x%jco(xnnz) = eqnnum+m
    END IF
    x%beg(eqnnum+1) = xnnz+1
    
!           Generate the rhs:
    
    RL(eqnnum) = Source*s
  END DO
  
END DO
x%nnz =xnnz

DEALLOCATE( h, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( xc, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( cw, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( ce, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( s1, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

RETURN
END SUBROUTINE arthura

END MODULE