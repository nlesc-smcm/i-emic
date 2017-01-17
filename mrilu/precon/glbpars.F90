MODULE m_glbpars

!
!-----------------------------------------------------------------------

!     Global parameters used in the computation of the
!     preconditioner and in the solution of the linear system.

!     The parameters get a default value during declaration.
!     The values of these
!     parameters can be changed by the SUBROUTINE 'iniglb'.


!     Description of the constant parameters:
!     =======================================
!     MachPrc   The machine precision, a small positive value.
!     Neglgbl   A small positive, negligible, value, smaller than the
!               machine precision 'MachPrec'.
!               Numbers with absolute value <= Neglgbl can be
!               considered zero.

DOUBLE PRECISION :: machprc, neglgbl
PARAMETER      ( machprc = 1.12D-16 , neglgbl = 1.0D-20 )


!     Description of the parameter in /glbpars/:
!     ==========================================
!     OutLev    Output level; the value of this parameter determines
!               the amount of output to the standard output file.
!               =  0  Only the error messages are written to standard
!                     output.
!               =  1  The output of level 0 plus:
!                     the warning messages.
!               =  2  The output of level 0 and 1 plus:
!                     the values of the input parameters.
!               =  3  The output of levels 0, 1 and 2 plus:
!                     the execution times.
!               =  4  The output of levels 0, 1, 2 and 3 plus:
!                     the sizes of the submatrices in the construction
!                     of the Multi Level Preconditioner, and
!                     the residuals during the iteration steps in the
!                     CG type solver.
!               >= 5  The output of levels 0, 1, 2, 3 and 4 plus:
!                     Extra output.

INTEGER :: outlev = 1

!-----------------------------------------------------------------------

END MODULE