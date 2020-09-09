      MODULE m_mat
        ! defines the location of the matrix
        ! replaces old common block file "mat.com" 

        ! originally in usr.com:
        real,    dimension(:,:,:,:,:,:), ALLOCATABLE :: Al

        ! originally in mat.com:
        real,    dimension(:), ALLOCATABLE :: coA
        integer, dimension(:), ALLOCATABLE :: jcoA 
        integer, dimension(:), ALLOCATABLE :: begA
        real,    dimension(:), ALLOCATABLE :: coB

        logical, dimension(:,:,:), ALLOCATABLE :: active
        ! used to keep track on active couplings

      END MODULE m_mat
