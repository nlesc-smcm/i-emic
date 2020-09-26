      MODULE m_bgsprec

!-----------------------------------------------------------------------
! Definition and description of block-Gauss-Seidel preconditioner
!

      USE m_build

!     Splits the matrix into the following structure
!     
!     [ Auv    0     Guv  BuvTS ]
!     [ 0      0     Gw   BwTS  ]
!     [ Duv    Dw    0    0     ]
!     [ BTSuv  BTSw  0    ATS   ]

!     

      TYPE bgsprec
         !  subblocks:
         TYPE (csrmatrix), POINTER :: Auv
         TYPE (csrmatrix), POINTER :: Guv
         TYPE (csrmatrix), POINTER :: BuvTS
         TYPE (csrmatrix), POINTER :: Gw
         TYPE (csrmatrix), POINTER :: BwTS
         TYPE (csrmatrix), POINTER :: Duv
         TYPE (csrmatrix), POINTER :: Dw
         TYPE (csrmatrix), POINTER :: BTSuv
         TYPE (csrmatrix), POINTER :: BTSw
         TYPE (csrmatrix), POINTER :: ATS
         !  dummy matrices
         TYPE (csrmatrix), POINTER :: Adumw
         TYPE (csrmatrix), POINTER :: Adump
         TYPE (csrmatrix), POINTER :: Bpdumw
         TYPE (csrmatrix), POINTER :: Bdumwp
         TYPE (csrmatrix), POINTER :: Buvdumw
         TYPE (csrmatrix), POINTER :: Bdumwuv
         !  depth-averaging  
         TYPE (csrmatrix), POINTER :: M1   ! chosen such that Gw*M1' = 0
!         TYPE (csrmatrix), POINTER :: M11  ! M1 = [M11 M12] 
!         TYPE (csrmatrix), POINTER :: M12  ! M12 square and invertible
         TYPE (csrmatrix), POINTER :: M2   ! chosen such that M2*Dw = 0
!         TYPE (csrmatrix), POINTER :: M21  ! M2 = [M21 M22]
!         TYPE (csrmatrix), POINTER :: M22  ! M22 square and invertible
         !  NB in case of stretched grid M2zp \= Mzp
         TYPE (csrmatrix), POINTER :: M1T   ! = M1'
!         TYPE (csrmatrix), POINTER :: M2T   ! = M2'
         TYPE (csrmatrix), POINTER :: GuvM1T  ! = Guv*M1'
         TYPE (csrmatrix), POINTER :: M2Duv   ! = M2*Duv
         !  in case of Modified Simpler approach
         TYPE (csrmatrix), POINTER :: Cp    ! = MDuv*inv(DMAuv)*MGuv
         !  in case of depth-averaged grad-div approach
         TYPE (csrmatrix), POINTER :: Mzuv
         TYPE (csrmatrix), POINTER :: MzuvT ! = Mzuv'
         TYPE (csrmatrix), POINTER :: MAuv  ! = Mzuv*Auv*Mzuv'
         TYPE (csrmatrix), POINTER :: MGuv  ! = Mzuv*Guv*Mzp'
         TYPE (csrmatrix), POINTER :: MDuv  ! = Mzp *Duv*Mzuv'
         !  grad-div stabilized matrix
         REAL                      :: omega
         TYPE (csrmatrix), POINTER :: MAuvGradDiv ! = MAuv + omega*MGuv*MDuv
         !  
         TYPE (csrmatrix), POINTER :: MCp   ! = MDuv*inv(DMAuv)*MGuv
         !  blocks of Ap
         TYPE (csrmatrix), POINTER :: Ap    ! square upper triangular part of Gw
         TYPE (csrmatrix), POINTER :: InvAp ! inverse of Ap (also upper tr.)
         TYPE (csrmatrix), POINTER :: Guv1  ! first d_w columns of Guv
         !  block of Aw
         TYPE (csrmatrix), POINTER :: Aw    ! square lower triangular part of Dw
         TYPE (csrmatrix), POINTER :: InvAw ! inverse of Aw
         TYPE (csrmatrix), POINTER :: Duv1  ! first d_w rows of Duv
         !  diagonal blocks of Auv
         TYPE (csrmatrix), POINTER :: DAuv    ! block-diagonal of Auv (blocksize=2)
         !  transformation of TS ! NOT used
         TYPE (csrmatrix), POINTER :: DTS        
         TYPE (csrmatrix), POINTER :: DTST        
         TYPE (csrmatrix), POINTER :: DTSinv
         TYPE (csrmatrix), POINTER :: DTSinvT
         !  preconditioners
         TYPE (prcmatrix), POINTER :: PAuv
         TYPE (csrmatrix), POINTER :: PDAuv
         TYPE (prcmatrix), POINTER :: PCp
         TYPE (prcmatrix), POINTER :: PATS
         !  preconditioners (2)
         TYPE (csrmatrix), POINTER :: PDMAuv
         TYPE (prcmatrix), POINTER :: PMAuv
         TYPE (prcmatrix), POINTER :: PMCp
         !  singular vectors of pressure
         REAL, DIMENSION(:), POINTER :: svp1
         REAL, DIMENSION(:), POINTER :: svp2
         REAL, DIMENSION(:), POINTER :: svp3
         REAL, DIMENSION(:), POINTER :: svp4
         !  transposed system or not ::
         LOGICAL                   :: transposed
      END TYPE bgsprec

      !  counter for spp iterations
      INTEGER :: sppiter
      INTEGER :: var = 0 ! value used and set in bgskit (applybilu)

      END MODULE m_bgsprec
