!#begindoc

MODULE m_inisol

CONTAINS

SUBROUTINE inisol (acgtype, amgmres, amaxnits, aredtol, aabstol)

USE m_solpars

INTEGER, INTENT(IN)                      :: acgtype
INTEGER, INTENT(IN)                      :: amgmres
INTEGER, INTENT(IN)                      :: amaxnits
DOUBLE PRECISION, INTENT(IN)             :: aredtol
DOUBLE PRECISION, INTENT(IN)             :: aabstol



!     INItialise parameters for the iterative SOLve routine.

!     Initialises the variables in solpars.
!     These variables are the parameters which control the solve
!     routine(s).

!     See also the description in the file  'solpars.F90'

!     Arguments:
!     ==========
!     aCGType  i   Type of Congruent Gradient solver:
!                  0=CG,
!                  1=BiCGSTAB,   using a left preconditioner,
!                  2=BiCGSTAB-R, using a Right preconditioner,
!                  3=GMRES(M)
!     aMgmres  i   M, the dimension of the subspace used in GMres.
!     aMaxNits i   Maximum number of iterations in CG-like solver.
!     aRedTol  i   Reduction Tolerance: The required factor to reduce
!                  residual of preconditioned system.
!                  For  bicgstabr, however, the required factor to
!                  reduce the residual of the original system.
!     aAbsTol  i   Absolute Tolerance:  Upper bound residual of the
!                  preconditioned system.
!                  For  bicgstabr, however, the upper bound of the
!                  relative residual of the original system.

!#enddoc

!     Global variables:
!     =================

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'inisol'

PRINT '(A, X, A)' , 'Entry:', rounam
#endif

cgtype  = acgtype
mgmres  = amgmres
maxnits = amaxnits
redtol  = aredtol
abstol  = aabstol

!     End of  inisol
END SUBROUTINE inisol

END MODULE