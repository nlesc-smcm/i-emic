!#begindoc
 
MODULE m_defvals

USE m_glbpars
USE m_prcpars
USE m_solpars
USE m_vispars

BLOCK DATA defvals

!     Define default values of the variables, used as parameters, in the
!     common blocks:
!        /glbpars/, /prcpars/, /solpars/ and /vispars/.
!     A description of these parameters can be found in the files:
!        glbpars.I90, prcpars.I90, solpars.I90 and vispars.I90

!#enddoc

!     Global Parameters:
!     ==================

!     Default values in /glbpars/:
DATA  outlev  /1/

!     Default values in /prcpars/:
DATA  cutmck  /.false./, scarow /.false./,  xactelm /.true./
DATA  gusfctr /1.00D0/
DATA  clsonce /.false./, nlsfctr /1.00D-1/
DATA  epsw    /1.00D-1/, elmfctr /2.00D-1/
DATA  redfctr /8.00D-1/
DATA  schtol  /0.0D0/
DATA  denslim /1.00D-2/, globfrac/0.0/,     locfrac /0.0/
DATA  sparslim /6.66D-1/
DATA  ilutype /9/,       droptol /1.00D-3/, cpivtol /8.75D-1/
DATA  compfctr /1.00D0/
DATA  lutol   /1.0D-10/, singlu /.false./

!     Default values in /solpars/:
DATA  cgtype  /1/,       maxnits /100/,     mgmres /10/
DATA  abstol  /1.0D-6/,  redtol /1.0D-6/

!     Default values in /vispars/:
DATA  nrowrtn /0/,       nscwrtn /0/
DATA  visnro  /0/,       visnsc /0/
DATA  visasc  /.false./, visildu /.false./, vislsc /.false./

END BLOCK DATA defvals

END MODULE
