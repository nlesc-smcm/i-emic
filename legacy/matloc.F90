	MODULE m_matloc
          ! define location of matrix, preconditioner and reordering
          USE m_build
          USE m_bgsprec
          TYPE (csrmatrix), POINTER :: ixA
          TYPE (bgsprec),   POINTER :: ixBGS
          TYPE (prcmatrix), POINTER :: ixPrc
          INTEGER,          DIMENSION(:),  ALLOCATABLE :: ord 
          REAL,             DIMENSION(:),  ALLOCATABLE :: svp1 
          REAL,             DIMENSION(:),  ALLOCATABLE :: svp2 
        END MODULE m_matloc
