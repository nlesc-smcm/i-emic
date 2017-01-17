!#begindoc
 
MODULE m_readpars

CONTAINS
 
SUBROUTINE readpars (readeq, blksiz, mint, probty, strtch, cvx, cvy, testx)

USE m_getunit
USE m_guerrmsg
USE m_rdfilenm
USE m_rdpblpars
USE m_rdprcpars
USE m_rdsolpars
USE m_rdvispars

INTEGER, INTENT(OUT)                  :: readeq
INTEGER, INTENT(OUT)                  :: blksiz
INTEGER, INTENT(OUT)                  :: mint
INTEGER, INTENT(OUT)                  :: probty
DOUBLE PRECISION, INTENT(OUT)         :: strtch
DOUBLE PRECISION, INTENT(OUT)         :: cvx
DOUBLE PRECISION, INTENT(OUT)         :: cvy
LOGICAL, INTENT(OUT)                  :: testx

!     Read the parameters for the multilevel preconditioner and for the
!     solver from a file.  The file name of this file is read from the
!     standard input file.

!     Arguments for the problem definition:
!     =====================================
!     ReadEq   o   Code number to Read or generate the linear system.
!     BlkSiz   o   Size of diagonal blocks in matrix linear system.
!     Mint     o   Number of intervals in X- and Y-direction.
!     ProbTy   o   Problem Type.
!     Strtch   o   Grid Stretching
!     CVX      o   Coefficient of  u_x
!     CVY      o   Coefficient of  u_y
!     TestX    o   Modify right-hand side for a given solution.

!     The values of the other parameters are directly stored into the
!     common blocks.  See also the subroutines:
!        iniprc, inisol  and  inivis

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'readpars'

INTEGER 		:: iunit, ier
CHARACTER (LEN=30) 	:: inpfnm


#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif


CALL getunit (iunit, ier)
IF ( ier /= 0 ) THEN
  CALL guerrmsg (rounam, ier)
  STOP 'program'
END IF

!     Read the name of the input file:
CALL rdfilenm ('Enter name input file:', inpfnm)

!     Open Input-file:

OPEN (iunit, FILE = inpfnm, STATUS = 'OLD')
REWIND (iunit)


!     Read the parameters which define the problem:

CALL rdpblpars (iunit, readeq, blksiz, mint, probty, strtch, cvx, cvy, testx)


!     Read the parameters to control the visualisation of the
!     (sub-)matrices during the construction of the Multi Level
!     Preconditioner, and initialise.

CALL rdvispars (iunit)


!     Read and check the parameters controlling
!        the preprocessing of the linear system, and
!        the construction of the Multi Level Preconditioner.

CALL rdprcpars (iunit)


!     Read the parameters to control the solution of the linear system
!     with a CG-like method, and initialise.

CALL rdsolpars (iunit)


CLOSE (iunit)

!     End of  readpars
END SUBROUTINE readpars

END MODULE