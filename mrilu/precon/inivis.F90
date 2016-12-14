!#begindoc
 
MODULE m_inivis

CONTAINS

SUBROUTINE inivis (avisasc, avisnsc, avislsc, avisildu, avisnro)

USE m_vispars

LOGICAL, INTENT(IN)                      :: avisasc
INTEGER, INTENT(IN)                      :: avisnsc
LOGICAL, INTENT(IN)                      :: avislsc
LOGICAL, INTENT(IN)                      :: avisildu
INTEGER, INTENT(IN)                      :: avisnro





!     INItialise parameters which influence the VISualisation of the
!     the matrices computed during the construction of the
!     preconditioner.

!     See also the description in the file  'vispars.F90'

!     Arguments:
!     ==========
!     aVisAsc  i   Visualise the original matrix A (possibly scaled).
!     aVisNSC  i   Visualise the first 'avisNSC' Schur-complement
!                  matrices.
!     aVisLSC  i   Visualise the last Schur-complement.
!     aVisILDU i   Visualize the Incomplete LDU factorization of the
!                  last Schur-complement.
!     aVisNRO  i   Visualise the first 'avisNRO' reordered matrices.

!#enddoc

!     Global variables:
!     =================

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'inivis'

PRINT '(A, X, A)' , 'Entry:', rounam
#endif

visasc   = avisasc
visnsc   = avisnsc
vislsc   = avislsc
visildu  = avisildu
visnro   = avisnro

nrowrtn  = 0
nscwrtn  = 0

!     End of  inivis
END SUBROUTINE inivis

END MODULE
