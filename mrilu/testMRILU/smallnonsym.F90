!=======================================================================

PROGRAM smallnonsym

!      4 Maart 1996.
!      4 Maart 2003 (Laatste aanpassing).
!      Testprogramma voor het oplossen van een stelsel lineaire
!      vergelijkingen
!         A*X=RL
!      Een meer algemene versie van dit programma staat in de
!      directory
!         ../testMRILU

!      Algemenere discretisatie op rooster met stretching
!      E.F.F.Botta@math.rug.nl

!      De grootte van de lokale arrays wordt gegeven door het aantal
!      roosterpunten 'Nm'. De waarde van deze variabele wordt gegeven
!      in een parameter statement in de routine 'eugen2d'.

!-----------------------------------------------------------------------

#ifndef WITH_UNION

#define prcmatrix  anymatrix
#define csrmatrix  anymatrix

#endif

USE m_build
USE m_matvec
USE m_iniglb
USE m_readpars
USE m_eugen2d

USE m_solprc
USE m_cmpprc


!     Local Parameters:
!     =================
CHARACTER (LEN=*),PARAMETER      :: rounam = 'smallnonsym'

INTEGER 					:: ier
INTEGER 					:: blksiz, n
INTEGER 					:: nblock, mint, probty
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE	:: dxx
DOUBLE PRECISION, DIMENSION(:), POINTER		:: RL, Xex
TYPE (csrmatrix), POINTER			:: x
LOGICAL 					:: testx
INTEGER 					:: readeq
INTEGER 					:: i
DOUBLE PRECISION 				:: strtch, cvx, cvy
DOUBLE PRECISION 				:: abserr
TYPE (prcmatrix), POINTER			:: Prc

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif


!     Initialize the global parameter 'OutLev':
CALL iniglb (4)


!     READ THE INPUT DATA and initialise the parameters for the
!     preconditioner:

CALL readpars (readeq, blksiz, mint, probty, strtch, cvx, cvy, testx)


!     Set up the matrix at [x], the right-hand
!     side at [RL], and the initial solution in [dxX]:

IF (readeq == 1) THEN
!        Determine the matrix and an exact solution from the convection-
!        diffusion equation on a 2D-grid:
  CALL eugen2d (probty, mint, strtch, cvx, cvy, n, x, Xex)
ELSE
  PRINT '(A, I11)' , 'Illegal code to initialise equation:', readeq
  STOP 'Error in  smallnonsym'
END IF

!     Set up the block sizes
nblock = n/blksiz

IF ( MOD(n,nblock) /= 0 ) THEN
  PRINT '(A, I10, 2X, A, I3)' , 'Number of equations', n,  &
      'NOT an integral multiple of block size', blksiz
END IF

PRINT '(/, 3(A, X, I10,/))' , 'Number of equations:', n,  &
    'Default block size: ', blksiz, 'Number of blocks:   ', nblock

ALLOCATE ( RL(1:n), STAT=ier)
IF (ier < 0)  GO TO 1000

!     Exact solution in 'Xex', determine right hand side in 'RL',
!     with   RL := A Xex:
RL = 0.0D0
CALL matvec (n, 1.0D0, csrtoany(x), Xex(1:n), RL)

!     Zero the initial solution:
dxx(1:x%n) = 0.0D0

!     Construct complete preconditioner and store location in 'Prc':
CALL cmpprc (blksiz, x, Prc)

!     Solve the linear system
!        A x = RL
!     using the complete preconditioner [Prc]:
CALL solprc (Prc, dxx, RL)

!     Test problem or problem derived from convection-diffusion
!     equation, with exact solution in Xex.

!     Calculate the absolute error.

abserr = MAXVAL( ABS((dxx(1:n)-Xex(1:n)) - (dxx(1))-Xex(1)) )

PRINT '(/, A, 1P, E22.16)', 'Absolute error:  ', abserr

!     Normal end:
STOP ': smallnonsym'

1000 CONTINUE
STOP ': smallnonsym'

!     End of  smallnonsym
END PROGRAM smallnonsym

