! North Atlantic
      integer, parameter :: n    =  80
      integer, parameter :: m    =  80
      integer, parameter :: l    =  16
      integer, parameter :: nun  =  6
      integer, parameter :: la   =  0
      integer, parameter :: ndim =  nun*n*m*(l+la)

      integer nid
      common /newton/ nid

      integer, parameter :: np   =  27	
      integer, parameter :: npar =  30	
      real         par
      common /par/ par(npar)

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
!     Par(20) determines the amount of 'pseudo' arc-length continuation.

      character*4 f_name ! run id.		! ATvS-VDA
      real dt_sc				! ATvS-VDA
      common /vda/ dt_sc, f_name		! ATvS-VDA
