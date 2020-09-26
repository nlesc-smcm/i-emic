	MODULE m_bgspars
          ! contains all parameters relevant for the solution process
          ! with the block-Gauss-Seidel preconditioner

          ! === PARAMETER FOR MRILU-OUTPUT CONTROL ===
          INTEGER, PARAMETER :: outlev = 0      ! value between 0 and 5 
          
          ! === PARAMETER FOR BGS/ILU ===
          LOGICAL, PARAMETER :: ilu   = .false. ! ILU parameter
          LOGICAL, PARAMETER :: xctsc = .false. ! compute "exact" schur complement?

          ! === PARAMETERS FOR THE KRYLOV METHODS ===
          ! extra: the solution parameters used in (f)gmres
          INTEGER, PARAMETER :: cgtype    = 3     ! select type of Krylov method 
          ! 1 : bicgstab (bad convergence)
          ! 2 : gmres    (acceptable performance)
          ! 3 : fgmres   (should be better)
          ! for the inner iterations always a variant of (f)gmres is used
          INTEGER            :: mgmres    = 30    ! number of iterations before restart
                                                  ! in (f)gmres
          INTEGER            :: maxnits   = 30    ! maximal number of iterations
                                                  ! should be a multiple of mgmres
          REAL               :: tolabs            ! set in bag.f or time.f
          REAL               :: tolred            ! set in bag.f or time.f
          ! nested iterations?
          LOGICAL, PARAMETER :: locgmres  =.false. ! how to solve Auv, ATS
                                                  ! .true. : with preconditioned gmres
                                                  ! .false.: apply preconditioner only once
          REAL,    PARAMETER :: loctolabs = 1E-6  ! tolabs for inner iterations
          REAL,    PARAMETER :: loctolred = 1E-6  ! tolred for inner iterations
          ! loctols for spp (always solved iteratively)
          REAL,    PARAMETER :: spploctolabs = 1E-6! tolabs for inner iterations spp
          REAL,    PARAMETER :: spploctolred = 1E-6! tolred for inner iterations spp
          ! to keep convergence the loctols for the spp have to be chosen small

          ! === PARAMETERS FOR THE SADDLE POINT PROBLEM ===
          CHARACTER (LEN=2),PARAMETER :: prectype  =  'AC'
          ! WS : Wathen/Sylvester preconditioner 
          ! ES : Elman/Sylvester preconditioner 
          ! GD : Grad-div stabilization preconditioner
          ! AC : Artificial compressibility preconditioner
          REAL,             PARAMETER :: sppomega  =  1.0D4
          ! for GD and AC: without scaling 1 appears to be appropriate
          !                with scaling it should be larger 100(?)
          !                depends strongly on the size of the problem
          ! for WS and ES: depends on scaling, but in general should be small
          !                with scaling 1D-2, without 1D-4/1D-5
          !                should not matter accoording to theory, nevertheless it does
          

          ! === MRILU PARAMETERS FOR MAuv(GradDiv) ===
          ! parameters for MAuv (if prectype = EW,ES) or
          ! its grad-div stabilized version MAuvGradDiv (if prectype = AC,GD)

          LOGICAL,          PARAMETER :: MAuv_cutmck   = .FALSE.
          LOGICAL,          PARAMETER :: MAuv_scarow   = .TRUE.
          LOGICAL,          PARAMETER :: MAuv_xactelm  = .TRUE. 
          LOGICAL,          PARAMETER :: MAuv_clsonce  = .FALSE.

          DOUBLE PRECISION, PARAMETER :: MAuv_nlsfctr  = 1.00D-1
          DOUBLE PRECISION            :: MAuv_epsw     = 1.00D-4
          DOUBLE PRECISION, PARAMETER :: MAuv_elmfctr  = 2.00D-1
          LOGICAL         , PARAMETER :: MAuv_gusmod   = .TRUE. 
          DOUBLE PRECISION, PARAMETER :: MAuv_gusfctr  = 1.00D00
          DOUBLE PRECISION, PARAMETER :: MAuv_redfctr  = 8.00D-1
          DOUBLE PRECISION, PARAMETER :: MAuv_schtol   = 0.0D0

          DOUBLE PRECISION, PARAMETER :: MAuv_denslim  = 1.00D-1
          DOUBLE PRECISION, PARAMETER :: MAuv_globfrac = 0.0
          DOUBLE PRECISION, PARAMETER :: MAuv_locfrac  = 0.1
          DOUBLE PRECISION, PARAMETER :: MAuv_sparslim = 6.66D-1
          
          INTEGER,          PARAMETER :: MAuv_ilutype  = 9
          DOUBLE PRECISION, PARAMETER :: MAuv_droptol  = 1.00D-7
          DOUBLE PRECISION, PARAMETER :: MAuv_compfct  = 1.00D0
          DOUBLE PRECISION, PARAMETER :: MAuv_cpivtol  = 8.75D-1
          DOUBLE PRECISION, PARAMETER :: MAuv_lutol    = 1.0D-10
          LOGICAL,          PARAMETER :: MAuv_singlu   = .TRUE.
      

          ! === MRILU PARAMETERS FOR Auv ===
          ! parameters for Auv
          LOGICAL,          PARAMETER :: Auv_cutmck   = .FALSE.
          LOGICAL,          PARAMETER :: Auv_scarow   = .TRUE.
          LOGICAL,          PARAMETER :: Auv_xactelm  = .TRUE. 
          LOGICAL,          PARAMETER :: Auv_clsonce  = .FALSE.

          DOUBLE PRECISION, PARAMETER :: Auv_nlsfctr  = 1.00D-1
          DOUBLE PRECISION            :: Auv_epsw     = 1.00D-4
          DOUBLE PRECISION, PARAMETER :: Auv_elmfctr  = 2.00D-1
          LOGICAL,          PARAMETER :: Auv_gusmod   = .TRUE.
          DOUBLE PRECISION, PARAMETER :: Auv_gusfctr  = 1.00D00
          DOUBLE PRECISION, PARAMETER :: Auv_redfctr  = 8.00D-1
          DOUBLE PRECISION, PARAMETER :: Auv_schtol   = 0.0D0
          
          DOUBLE PRECISION, PARAMETER :: Auv_denslim  = 3.00D-1
          DOUBLE PRECISION, PARAMETER :: Auv_globfrac = 0.0
          DOUBLE PRECISION, PARAMETER :: Auv_locfrac  = 0.1
          DOUBLE PRECISION, PARAMETER :: Auv_sparslim = 6.66D-1
          
          INTEGER,          PARAMETER :: Auv_ilutype  = 9
          DOUBLE PRECISION, PARAMETER :: Auv_droptol  = 1.00D-4
          DOUBLE PRECISION, PARAMETER :: Auv_compfct  = 1.00D0
          DOUBLE PRECISION, PARAMETER :: Auv_cpivtol  = 8.75D-1
          DOUBLE PRECISION, PARAMETER :: Auv_lutol    = 1.0D-10
          LOGICAL,          PARAMETER :: Auv_singlu   = .FALSE.


          ! === MRILU PARAMETER FOR ATS ===
          ! parameters for ATS
          LOGICAL,          PARAMETER :: ATS_cutmck   = .FALSE.
          LOGICAL,          PARAMETER :: ATS_scarow   = .TRUE.
          LOGICAL,          PARAMETER :: ATS_xactelm  = .TRUE.
          LOGICAL,          PARAMETER :: ATS_clsonce  = .TRUE.

          DOUBLE PRECISION, PARAMETER :: ATS_nlsfctr  = 1.00D-1
          DOUBLE PRECISION            :: ATS_epsw     = 3.00D-2
          DOUBLE PRECISION, PARAMETER :: ATS_elmfctr  = 2.00D-1
          LOGICAL,          PARAMETER :: ATS_gusmod   = .TRUE.
          DOUBLE PRECISION, PARAMETER :: ATS_gusfctr  = 0.95D0
          DOUBLE PRECISION, PARAMETER :: ATS_redfctr  = 1.00D0
          DOUBLE PRECISION, PARAMETER :: ATS_schtol   = 0.0D0
          
          DOUBLE PRECISION, PARAMETER :: ATS_denslim  = 9.00D-1
          DOUBLE PRECISION, PARAMETER :: ATS_globfrac = 0.0
          DOUBLE PRECISION, PARAMETER :: ATS_locfrac  = 0.01
          DOUBLE PRECISION, PARAMETER :: ATS_sparslim = 9.50D-1
          
          INTEGER,          PARAMETER :: ATS_ilutype  = 9
          DOUBLE PRECISION, PARAMETER :: ATS_droptol  = 3.00D-5
          DOUBLE PRECISION, PARAMETER :: ATS_compfct  = 1.00D0
          DOUBLE PRECISION, PARAMETER :: ATS_cpivtol  = 8.75D-1
          DOUBLE PRECISION, PARAMETER :: ATS_lutol    = 1.0D-10
          LOGICAL,          PARAMETER :: ATS_singlu   = .TRUE.

        END MODULE m_bgspars
