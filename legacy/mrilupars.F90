      MODULE m_mrilupars

        ! prcpars
        REAL    :: compfctr = 1.0E0
        REAL    :: cpivtol  = 1.0E-1
        REAL    :: denslim  = 1.0E-3
        REAL    :: droptol  = 5.0E-4
        REAL    :: epsw     = 1.0E-7
        REAL    :: elmfctr  = 0.0E0
        REAL    :: globfrac = 1.0
        LOGICAL :: gusmod   = .false.
        REAL    :: gusfctr  = 1.0
        REAL    :: locfrac  = 1.0
        REAL    :: lutol    = 0.0E-12
        REAL    :: nlsfctr  = 1.0E0
        REAL    :: redfctr  = 0.8E0
        REAL    :: schtol   = 1.0E-6
        REAL    :: sparslim = 6.75E-1
  
        INTEGER :: ilutype = 9
        
        LOGICAL :: clsonce = .false.
        LOGICAL :: cutmck  = .false.
        LOGICAL :: scarow  = .true.
        LOGICAL :: singlU  = .true.
        LOGICAL :: xactelm = .true.
        
        ! solpars
        INTEGER :: cgtype1 = 3
        INTEGER :: mgmres1 = 60
        INTEGER :: maxnits1= 120

        REAL    :: redtol1 = 1.0E-7
        REAL    :: abstol1 = 1.0E-6

      END MODULE m_mrilupars
      
