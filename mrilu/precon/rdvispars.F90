!#begindoc
 
MODULE m_rdvispars

CONTAINS

SUBROUTINE rdvispars (iunit)

USE m_inivis
USE m_dump

INTEGER, INTENT(IN)                  :: iunit

!     Read the parameters to control the visualisation of the
!     (sub-)matrices during the construction of the Multi Level
!     Preconditioner, and initialise common block
!     /vispars/.

!     The expected format of the input file consists of the lines with:
!        A dummy line
!        Letter (T/F)          VisAsc.  Visualise the (scaled) matrix A.
!        Integer (0,1,2,...)   VisNSC.  Visualise the first 'visNSC'
!                              Schur-complement matrices.
!        Letter (T/F)          VisLSC.  Visualise the last Schur-
!                              complement.
!        Letter (T/F)          VisILDU.  Visualize the Incomplete LDU
!                              factorization of the last Schur-
!                              complement.
!        Integer (0,1,2,...)   VisNRO   Visualise the first 'visNRO'
!                              ReOrdered matrices.

!     Argument:
!     =========
!     iunit    i   Unit number of the open input file containing the
!                  input parameters.

!#enddoc

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'rdvispars'
#endif

CHARACTER (LEN=80) :: line

INTEGER :: visnro, visnsc
LOGICAL :: visasc, vislsc, visildu

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Skip the first line:
READ (iunit, '(A80)') line

!     Visualisation parameters:
READ (iunit, *) visasc

READ (iunit, *) visnsc

IF (visnsc < 0) CALL dump(__FILE__,__LINE__, 'Negative value VisNSC NOT allowed!')

READ (iunit, *) vislsc
READ (iunit, *) visildu

READ (iunit, *) visnro

IF (visnro < 0)  CALL dump(__FILE__,__LINE__, 'Negative value VisNRO NOT allowed!')

CALL inivis (visasc, visnsc, vislsc, visildu, visnro)


!     Normal return:
RETURN


!     Fatal error:
1000 CONTINUE
STOP 'in rdvispars: Wrong input value!'

!     End of  rdvispars
END SUBROUTINE rdvispars

END MODULE