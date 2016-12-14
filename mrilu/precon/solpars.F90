 
!-----------------------------------------------------------------------

!     Solver parameters.

!     The parameters get a default value
!     during declaration. The values of these
!     parameters can be changed by the SUBROUTINE 'inisol'.


!     Description of the parameters in /solpars/:
!     ===========================================
!     CGType    Type of Congruent Gradient solver:
!               0=CG,
!               1=BiCGSTAB,   using a left preconditioner,
!               2=BiCGSTAB-R, using a Right preconditioner,
!               3=GMRES(M)
!     Mgmres    M, the dimension of the subspace used in GMres.
!     MaxNIts   Maximum number of iterations in CG-like solver.
!     RedTol    Reduction Tolerance: The required factor to reduce
!               residual of preconditioned system.
!               For  bicgstabr, however, the required factor to reduce
!               residual of the original system.
!     AbsTol    Absolute Tolerance:  Upper bound residual of the
!               preconditioned system.
!               For  bicgstabr, however, the upper bound of the
!               relative residual of the original system.

MODULE m_solpars

DOUBLE PRECISION :: 		&
	redtol	=1.0D-6		,&
	abstol	=1.0D-6

INTEGER ::			&
	cgtype	=1		,&
	mgmres	=10		,&
	maxnits	=100

!COMMON  /solpars/ redtol, abstol, cgtype, mgmres, maxnits

!SAVE /solpars/

END MODULE

!-----------------------------------------------------------------------

