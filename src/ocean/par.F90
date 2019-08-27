!! A Fortran 90 replacement for the file 'par.com'.
!! simply replace "include 'par.com'" by "use m_par"
!! and everything should be fine. Note that the dimensions
!! n,m,l are now set in the function usr::init instead of
!! being defined here.
module m_par
  ! user parameters (previously located in usr.com)...........................
#ifdef DATA_DIR
  character(len=*), parameter ::  topdir = DATA_DIR
  character(len=*), parameter ::  rundir = ''
#else
#error no guess for data directory available!
#endif
  real,    parameter :: pi    =  3.14159265358979323846
  integer, parameter :: SLIP =  -1 ! noslip boundary SLIP = -1
  !! previously par.com......................................................
  !  number of unknowns
  integer, parameter :: nun  =  6

  integer :: nid

  ! size of stencil/neighbourhood:
  ! +----------++-------++----------+
  ! | 12 15 18 || 3 6 9 || 21 24 27 |
  ! | 11 14 17 || 2 5 8 || 20 23 26 |
  ! | 10 13 16 || 1 4 7 || 19 22 25 |
  ! |  below   || center||  above   |
  ! +----------++-------++----------+
  integer, parameter :: np   =  27

  ! number of adjustable parameters:
  integer, parameter :: npar =  30

  ! array containing all these parameters:
  real, dimension(npar) :: par

  ! enumeration to identify the parameters:
  integer, parameter :: AL_T   =  1
  integer, parameter :: RAYL   =  2
  integer, parameter :: EK_V   =  3
  integer, parameter :: EK_H   =  4
  integer, parameter :: ROSB   =  5
  integer, parameter :: MIXP   =  6
  integer, parameter :: RESC   =  7
  integer, parameter :: SPL1   =  8
  integer, parameter :: HMTP   =  9
  integer, parameter :: SUNP   = 10
  integer, parameter :: PE_H   = 11
  integer, parameter :: PE_V   = 12
  integer, parameter :: P_VC   = 13
  integer, parameter :: LAMB   = 14
  integer, parameter :: SALT   = 15
  integer, parameter :: WIND   = 16
  integer, parameter :: TEMP   = 17
  integer, parameter :: BIOT   = 18
  integer, parameter :: COMB   = 19
  integer, parameter :: ARCL   = 20
  integer, parameter :: NLES   = 21
  integer, parameter :: IFRICB = 22
  integer, parameter :: CONT   = 23
  integer, parameter :: ENER   = 24
  integer, parameter :: ALPC   = 25
  integer, parameter :: CMPR   = 26
  integer, parameter :: FPER   = 27
  integer, parameter :: SPER   = 28
  integer, parameter :: MKAP   = 29
  integer, parameter :: SPL2   = 30

  ! enumeration to identify the unknowns
  integer, parameter :: UU = 1
  integer, parameter :: VV = 2
  integer, parameter :: WW = 3
  integer, parameter :: PP = 4
  integer, parameter :: TT = 5
  integer, parameter :: SS = 6

  ! enumeration to identify grid point type
  integer, parameter :: OCEAN = 0
  integer, parameter :: LAND  = 1
  integer, parameter :: WATER = 2
  integer, parameter :: PERIO = 3

end module m_par
