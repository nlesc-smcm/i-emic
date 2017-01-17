!#begindoc

MODULE m_iniglb 

CONTAINS

SUBROUTINE iniglb (aoutlev)

USE m_glbpars

INTEGER, INTENT(IN)                      :: aoutlev


!     INItialise GLoBal parameters.


!     See also the description in the file  'glbpars.F90'

!     Arguments:
!     ==========
!     aOutLev  i   Output level; the value of this parameter determines
!                  the amount of output to the standard output file.
!                  =  0  Only the error messages are written to standard
!                        output.
!                  =  1  The output of level 0 plus:
!                        the warning messages.
!                  =  2  The output of level 0 and 1 plus:
!                        the values of the input parameters.
!                  =  3  The output of levels 0, 1 and 2 plus:
!                        the execution times.
!                  =  4  The output of levels 0, 1, 2 and 3 plus:
!                        the sizes of the submatrices in the construction
!                        of the Multi Level Preconditioner, and
!                        the residuals during the iteration steps in the
!                        CG type solver.
!                  >= 5  The output of levels 0, 1, 2, 3 and 4 plus:
!                        Extra output.

!#enddoc

!     Global variables:
!     =================

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'iniglb'

PRINT '(A, X, A)' , 'Entry:', rounam
#endif

outlev = aoutlev


!     End of  iniglb
END SUBROUTINE iniglb

END MODULE