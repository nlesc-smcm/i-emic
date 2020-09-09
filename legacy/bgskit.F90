      MODULE m_bgskit
!
      CONTAINS

!**************************************************************
      SUBROUTINE build_BGSprec(ixA,ixBGS,ord,svp1,svp2,type)
!     Builds the block-Gauss-Seidel preconditioner
      USE m_build
      USE m_bgsprec
      USE m_sparsekit
      USE m_bgspars
      USE m_wfree
      IMPLICIT none
!      INCLUDE 'par.com'
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER          :: ixA
      TYPE (bgsprec), POINTER            :: ixBGS
      INTEGER, DIMENSION(:)              :: ord
      REAL, DIMENSION(:)                 :: svp1,svp2
      INTEGER :: type ! 1 = depth averaged spp, 2 = modified simpler
!     LOCAL
      INTEGER, DIMENSION(:), ALLOCATABLE :: iord
      INTEGER                            :: dumw, dump, ier
!     construct subblocks and depth-averaged system
      WRITE (6,*) "Constructing block-Incomplete LU prec for Jacobian"
      ! 1) construct block-ordering and reoderder matrix
      ! 1a) construct general block ordering
      IF (outlev.ge.1) WRITE (6,*) "Constructing block ordering"
!      CALL write_sparse_matrix(ixA,"A")
      CALL block_ordering(ord)
      ! 1b) detect dummy w and p points
      IF (outlev.ge.1) WRITE (6,*) "Detecting dummy w and p points"
      CALL detect_dummies(ixA,ord,dumw,dump)
      ! 1c) construct inverse ordering
      IF (outlev.ge.1) WRITE (6,*) "Constructing inverse ordering"
      ALLOCATE(iord(ixA%n))
      CALL invert_ordering(ord,iord)
      ! 1d) reorder matrix
      IF (outlev.ge.1) WRITE (6,*) "Reordering matrix"
      CALL reorder_matrix(ixA,ord,iord)
      ! 2) construct building blocks of the preconditioner
      ! 2a) split the matrix
      IF (outlev.ge.1) WRITE (6,*) "Splitting matrix"
      CALL split_matrix(ixA,ixBGS,dumw,dump)
!      CALL transformTS(ixBGS)
      IF (type.eq.1) THEN
         ! 2c) compute depth-averaging matrices
         IF (outlev.ge.1) WRITE (6,*) "Computing depth-averaged equations"
         CALL depth_averaging(ixBGS,type)
         ! 2d) grad-div stabilization
         IF (outlev.ge.1) WRITE (6,*) "Computing Grad-Div stabilized matrix"
         CALL grad_div_stabilization(ixBGS,sppomega)
      ELSEIF (type.eq.2) THEN
         ! 2c) compute depth-averaging matrices
         IF (outlev.ge.1) WRITE (6,*) "Computing depth-averaged equations"
         CALL depth_averaging(ixBGS,type)
         ! 2d) modified simpler preconditioner
         IF (outlev.ge.1) WRITE (6,*) "Computing Modified Simpler preconditioner"
         CALL compute_simpler(ixBGS)
      ELSE
         STOP "Unknown type of preconditioner, change prec in usr.com"
      END IF
      ! 2e) build diagonal blocks Ap and Aw
      IF (outlev.ge.1) WRITE (6,*) "Building diagonal blocks Ap and Aw"
      CALL build_diag_blocks(ixBGS)
      ! 2f) construct singular vectors of the pressure
      IF (outlev.ge.1) WRITE (6,*) "Constructing singular vectors for the pressure"
      IF (ixBGS%transposed) THEN
         CALL build_svp3(ixBGS)
      ELSE
         CALL build_svp(ixBGS)
      END IF      
!      CALL frs_svps(ixBGS,1)
!      CALL build_svp(ixBGS)
!      CALL build_svp2(ixBGS)
      IF (outlev.ge.1) WRITE (6,*) "Computing MRILU factorizations"
      ! 3) compute LU factorizations for Auv, ATS and MAuvGradDiv
      CALL compute_precons(ixBGS,type)
      ! 4) restore original ordering of the matrix
      IF (outlev.ge.1) WRITE (6,*) "Reorder matrix back"
      CALL reorder_matrix(ixA,iord,ord)
      CALL compute_big_svps(ixBGS,svp1,svp2,dump)
      svp1 = svp1(iord)
      svp2 = svp2(iord)      
      DEALLOCATE(iord)
      WRITE (6,*) "Construction of preconditioner finished"

      END SUBROUTINE build_BGSprec

!**************************************************************
      SUBROUTINE block_ordering(ord)
!     Builds the blockordering (uv,w,p,T,S) and its inverse
      IMPLICIT none
      INCLUDE 'par.com'
!     IMPORT/EXPORT
      INTEGER, DIMENSION(:) :: ord
!     LOCAL
      integer NN,i
!
      DO i = 1,ndim
         ord(i) = i
      END DO
      END SUBROUTINE block_ordering

!**************************************************************
      SUBROUTINE detect_dummies(ixA,ord,dumw,dump)
!     Changes the ordering 'ord' such that the dummy w's and p's 
!     and the free-surface w's are in front
      USE m_build
      USE m_bgspars
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER                  :: ixA
      INTEGER, DIMENSION(:),INTENT(INOUT),TARGET :: ord
      INTEGER                                    :: dumw, dump
!     LOCAL
      INTEGER i,ier,dim,nr,nnzA,j,row,k
      INTEGER, DIMENSION(:), ALLOCATABLE  :: subord
      INTEGER, DIMENSION(:), POINTER      :: pointord
      LOGICAL, DIMENSION(:), ALLOCATABLE  :: flag
!     
      dim = ixA%n
      nr = dim/6
      j = 0
      ALLOCATE(subord(nr),flag(nr))
      flag = .false.
      ! detect dummy w's
      ! Dummy w's and p's can be detected because they have a
      ! nonzero entry on the diagonal, normal w's and p's don't.
      DO i = 1,nr
         row = ord(2*nr+i)
         loop1: DO k = ixA%beg(row),ixA%beg(row+1)-1
            IF ((ixA%jco(k)).EQ.row) THEN
               j = j+1
               subord(j) = i
               flag(i) = .true.
               EXIT loop1
            END IF
         END DO loop1
      END DO
      dumw = j
      IF (outlev.ge.1) WRITE(6,"("" Number of dummy ws:"",I7)") dumw      
      DO i = 1,nr
         IF (.NOT.flag(i)) THEN
            j = j+1
            subord(j) = i
         END IF
      END DO
      IF (j.NE.nr) THEN
         IF (outlev.ge.1) WRITE(6,"("" Empty w-rows detected!:"",I7)") nr-j
      END IF
      pointord => ord(2*nr+1:3*nr)
      pointord = pointord(subord)

      ! detect dummie p's
      j = 0
      flag = .false.
      DO i = 1,nr
         row = ord(3*nr+i)
         if (ixA%beg(row).eq.ixA%beg(row+1)) then
            WRITE(*,*) "WARNING empty p-row detected"
         END if
         loop: DO k = ixA%beg(row),ixA%beg(row+1)-1
            IF ((ixA%jco(k)).EQ.row) THEN
               j = j+1
               subord(j) = i
               flag(i) = .true.
               EXIT loop
            END IF
         END DO loop
      END DO
      dump = j
      IF (outlev.ge.1) WRITE(6,"("" Number of dummy ps:"",I7)") dump
      DO i = 1,nr
         IF (.NOT.flag(i)) THEN
            j = j+1
            subord(j) = i
         END IF
      END DO
      IF (j.NE.nr) THEN
         IF (outlev.ge.1) WRITE(6,"("" Empty p-rows detected!:"",I7)") nr-j
      END IF
      pointord => ord(3*nr+1:4*nr)
      pointord = pointord(subord)
      DEALLOCATE(subord)
      END SUBROUTINE detect_dummies

!**************************************************************
      SUBROUTINE split_matrix(ixA,ixBGS,dumw,dump)
!     Stores the following blocks of the matrix in 
!     the block-Gauss-Seidel preconditioner
!     
!     [ Auv    0       0     0      Guv     BuvTS ]
!     [ 0      Adumw   0     0      Bdumwp  0     ]
!     [ 0      0       0     0      Gw      BwTS  ]
!     [ 0      0       0     Adump  0       0     ]
!     [ Duv    Bpdumw  Dw    0      0       0     ]
!     [ BTSuv  0       BTSw  0      0       ATS   ]
!     
      USE m_build
      USE m_bgsprec
      USE m_bgspars
      USE m_sparsekit
      USE m_wfree
      IMPLICIT none
      INCLUDE 'usr.com'
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER    :: ixA,ixB
      TYPE (bgsprec), POINTER      :: ixBGS
      INTEGER                      :: dumw,dump
!     LOCAL
      INTEGER :: nr, nnzA,nA,ier

      ! allocate preconditioner
      ALLOCATE(ixBGS, STAT=ier)
      nr = n*m*(l+la)
      nA = ixA%n
      ! extract the convection-diffusion-Coriolis equation
      CALL extract_submatrix(ixA,ixBGS%Auv  ,1          ,2*nr,1          ,2*nr)
      CALL extract_submatrix(ixA,ixBGS%Guv  ,1          ,2*nr,3*nr+dump+1,4*nr)
      CALL extract_submatrix(ixA,ixBGS%BuvTS,1          ,2*nr,4*nr+1,     6*nr)
      ! extract the hydrostatic pressure equation
      CALL extract_submatrix(ixA,ixBGS%Gw   ,2*nr+dumw+1,3*nr,3*nr+dump+1,4*nr)
      CALL extract_submatrix(ixA,ixBGS%BwTS ,2*nr+dumw+1,3*nr,4*nr+1     ,6*nr)
      ! extract the continuity equation
      CALL extract_submatrix(ixA,ixBGS%Duv  ,3*nr+dump+1,4*nr,1     ,2*nr)
      CALL extract_submatrix(ixA,ixBGS%Dw   ,3*nr+dump+1,4*nr,2*nr+dumw+1,3*nr)
      ! multiply continuity equation with -1
      ixBGS%Duv%co = - ixBGS%Duv%co
      ixBGS%Dw %co = - ixBGS%Dw %co
      ! extract the convection-diffusion equation of heat and salt
      CALL extract_submatrix(ixA,ixBGS%BTSuv,4*nr+1,6*nr,1     ,2*nr)
      CALL extract_submatrix(ixA,ixBGS%BTSw ,4*nr+1,6*nr,2*nr+dumw+1,3*nr)
      CALL extract_submatrix(ixA,ixBGS%ATS  ,4*nr+1,6*nr,4*nr+1,6*nr)
      ! extract the dummy blocks
      CALL extract_submatrix(ixA,ixBGS%Adumw,2*nr+1,2*nr+dumw,2*nr+1,2*nr+dumw)
      CALL extract_submatrix(ixA,ixBGS%Adump,3*nr+1,3*nr+dump,3*nr+1,3*nr+dump)
      CALL extract_submatrix(ixA,ixBGS%Bpdumw,3*nr+dump+1,4*nr,2*nr+1,2*nr+dumw)
      CALL extract_submatrix(ixA,ixBGS%Bdumwp,2*nr+1,2*nr+dumw,3*nr+dump+1,4*nr)
      CALL extract_submatrix(ixA,ixBGS%Buvdumw,1,2*nr,2*nr+1,2*nr+dumw)
      CALL extract_submatrix(ixA,ixBGS%Bdumwuv,2*nr+1,2*nr+dumw,1,2*nr)
      ! multiply continuity equation with -1
      ixBGS%Bpdumw%co = - ixBGS%Bpdumw%co
      !zero check?
      nnzA = 0
      nnzA = nnzA + ixBGS%Auv%nnz    
      nnzA = nnzA + ixBGS%Guv%nnz    
      nnzA = nnzA + ixBGS%BuvTS%nnz    
      nnzA = nnzA + ixBGS%Gw%nnz    
      nnzA = nnzA + ixBGS%BwTS%nnz

      nnzA = nnzA + ixBGS%Duv%nnz   
      nnzA = nnzA + ixBGS%Dw%nnz

      nnzA = nnzA + ixBGS%BTSuv%nnz   
      nnzA = nnzA + ixBGS%BTSw%nnz    
      nnzA = nnzA + ixBGS%ATS%nnz     

      nnzA = nnzA + ixBGS%Adumw%nnz  
      nnzA = nnzA + ixBGS%Adump%nnz  
      nnzA = nnzA + ixBGS%Bpdumw%nnz 
      nnzA = nnzA + ixBGS%Bdumwp%nnz 
      
      IF (frs) THEN
         IF (ixBGS%Bpdumw%nnz.LT.ixBGS%Bdumwp%nnz) THEN
            ixBGS%transposed = .false.
         ELSE
            ixBGS%transposed = .true.
         END IF
      ELSE
         IF (ixBGS%Bpdumw%nnz.GT.0) THEN
            ixBGS%transposed = .false.
         ELSEIF (ixBGS%Bdumwp%nnz.GT.0) THEN
            ixBGS%transposed = .true.
         ELSE
            ixBGS%transposed = .false.
         END IF
      END IF
      IF ((outlev.ge.1).AND.ixBGS%transposed) THEN
         WRITE(*,*) "Solving transposed system"
      END IF

      nnzA = ixA%nnz - nnzA
      IF (nnzA.NE.0) THEN
         CALL unexpected_nnz(ixA,nnzA,dumw,dump,ixBGS%transposed)
      END IF

      END SUBROUTINE split_matrix

!**************************************************************
      SUBROUTINE unexpected_nnz(ixA,nnzA,dumw,dump,trans)
!     research and warning in case of unexpected nonzeros
      USE m_build
      USE m_bgsprec
      USE m_sparsekit
      USE m_bgspars
      USE m_wfree
      IMPLICIT none
      INCLUDE 'par.com'
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER    :: ixA
      INTEGER                      :: nnzA, dumw, dump
      LOGICAL                      :: trans
!     LOCAL
      TYPE (csrmatrix), POINTER    :: ixB
      INTEGER                      :: nr
      
      IF (trans) THEN
         nr = n*m*(l+la)
         CALL extract_submatrix(ixA,ixB,2*nr+dumw+1,3*nr,1,2*nr+1)
	 IF (outlev.gt.1) WRITE(6,*) "real w-coupling to uv's  :", ixB%nnz
         nnzA = nnzA - ixB%nnz
         CALL csrfree(ixB)
         CALL extract_submatrix(ixA,ixB,2*nr+1,2*nr+dumw,1,2*nr)
         IF (outlev.gt.1) WRITE(6,*) "dummy w coupling to uv   :", ixB%nnz
         nnzA = nnzA - ixB%nnz
         IF (outlev.gt.1) WRITE(6,*) "left over unexp.nonzeros :", nnzA
         CALL csrfree(ixB)
      ELSE
         nr = n*m*(l+la)
         CALL extract_submatrix(ixA,ixB,1,2*nr+1,2*nr+dumw+1,3*nr)
         IF (outlev.gt.1) WRITE(6,*) "uv-coupling to real w's  :", ixB%nnz
         nnzA = nnzA - ixB%nnz
         CALL csrfree(ixB)         
         CALL extract_submatrix(ixA,ixB,1,2*nr,2*nr+1,2*nr+dumw)
         IF (outlev.gt.1) WRITE(6,*) "uv-coupling to dummy w's :", ixB%nnz
         nnzA = nnzA - ixB%nnz
         IF (outlev.gt.1) WRITE(6,*) "left over unexp.nonzeros :", nnzA
         CALL csrfree(ixB)
      END IF
      
      IF (nnzA > 0) THEN
         WRITE(99,*) " "
         WRITE(99,*) "WARNING: Matrix contains a number of unexpected nonzeros: ", nnzA
         WRITE(99,*) "You probably changed the equations, the preconditioner might fail."
         WRITE(99,*) " "
      END IF
      END SUBROUTINE unexpected_nnz

!**************************************************************
      SUBROUTINE depth_averaging(ixBGS,type)
!     constructs the depth-averaging blocks of the preconditioner
!     and computes the depth-average saddle-point equation
      USE m_bgsprec
      USE m_bgspars
      USE m_sparsekit
      ! mrilu modules:
      USE m_build
      USE m_wfree
      USE m_wacsr
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (bgsprec),   POINTER :: ixBGS
      INTEGER :: type ! 1 = depth averaged spp, 2 = modified simpler
!     LOCAL
      TYPE (csrmatrix), POINTER :: ixC,ixB,ixDwT
      INTEGER                   :: ier, i,j,k
      REAL                      :: nrm
!
      ! copmute Mzp: the singular vectors of Gw
      CALL build_singular_matrix(ixBGS%Gw,ixBGS%Gw%n,ixBGS%Dw%n,ixBGS%M1)
      ! compute M2zp: the singular vectors of Dw'
      CALL transpose_csr(ixBGS%Dw,ixBGS%Dw%n,ixBGS%Gw%n,ixDwT)
      CALL build_singular_matrix(ixDwT,ixBGS%Gw%n,ixBGS%Dw%n,ixBGS%M2)
      CALL csrfree(ixDwT)
      ! compute Mzp'
      CALL transpose_csr(ixBGS%M1,ixBGS%M1%n,ixBGS%Dw%n,ixBGS%M1T)
!      CALL transpose_csr(ixBGS%M2,ixBGS%M2%n,ixBGS%Gw%n,ixBGS%M2T)

      ! check if Gw*M1' = 0
      CALL mult(ixBGS%Gw,ixBGS%M1T, ixC,2*ixBGS%M1T%nnz)
      nrm = dsqrt(dot_product(ixC%co,ixC%co))
      IF ((nrm.GT.1D-12).AND.(outlev.ge.1)) THEN 
         WRITE(6,"("" in depth_averaging: norm(Gw*Mzp') = "",D16.6)") nrm
      END IF
      CALL csrfree(ixC)

      ! check if M2*Dw = 0
      CALL mult(ixBGS%M2,ixBGS%Dw, ixC,2*ixBGS%M2%nnz)
      nrm = dsqrt(dot_product(ixC%co,ixC%co))
      IF ((nrm.GT.1D-12).AND.(outlev.ge.1)) THEN 
         WRITE(6,"("" in depth_averaging: norm(M2zp*Dw) = "",D16.6)") nrm
      END IF
      CALL csrfree(ixC)

      ! construct M2Duv and GuvM1T
      CALL mult(ixBGS%Guv, ixBGS%M1T, ixBGS%GuvM1T, 2*ixBGS%Guv%nnz)
      CALL mult(ixBGS%M2,  ixBGS%Duv, ixBGS%M2Duv,  2*ixBGS%Duv%nnz)
      IF (type.eq.1) THEN
         CALL build_Mzuv(ixBGS)
         CALL transpose_csr(ixBGS%Mzuv, ixBGS%Mzuv%n,ixBGS%Auv%n,ixBGS%MzuvT)
         ! finally construct the depth average saddle point equation
         CALL depth_average_block(ixBGS%Auv,ixBGS%Mzuv,ixBGS%MzuvT,ixBGS%MAuv)
         CALL mult(ixBGS%Mzuv,  ixBGS%GuvM1T, ixBGS%MGuv, 2*ixBGS%GuvM1T%nnz)
         CALL mult(ixBGS%M2Duv, ixBGS%MzuvT,  ixBGS%MDuv, 2*ixBGS%M2Duv%nnz)
      END IF
      END SUBROUTINE depth_averaging

!**************************************************************
      SUBROUTINE build_Mzuv(ixBGS)
!     constructs the depth-averaging matrices Mzp and Mzuv
!     such that the matrices are normalized : Mzp*Mzp'=Id, Mzuv*Mzuv'=Id
      USE m_build
      USE m_bgsprec
      USE m_wacsr
!      USE m_wecsr
      USE m_wfree
      USE m_sparsekit
      IMPLICIT none
      INCLUDE 'usr.com'
!     IMPORT/EXPORT
      TYPE (bgsprec),        POINTER :: ixBGS
!     LOCAL
      TYPE (csrmatrix),      POINTER :: ix1, ix2
      INTEGER                        :: nr,nl,ier,i,j,k, dumuv,svs
      INTEGER                        :: rowu, rowv, colu, colv, idx
!     

      nr = ixBGS%Auv%n
      nr = nr   ! real dimension of the problem
      nl = nr/(l+la) ! dimension of one layer
      CALL wacsr(nl,nr,ixBGS%Mzuv)
      ixBGS%Mzuv%beg(nl+1) = nr+1
      DO i = 1,nl
         ixBGS%Mzuv%beg(i) = (i-1)*l + 1
         DO k = 1,(l)
            idx = (i-1)*l+k
            ixBGS%Mzuv%jco(idx) = i + (k-1)*nl 
         END DO
      END DO
      
!      ixBGS%Mzuv%co = 1
      ixBGS%Mzuv%co = 1D0/sqrt(real(l))
      ! set coefficients of landpoints to 0
      DO i = 1,n-1
      DO j = 1,m-1
      loop:DO k = 1,l-1
         IF ((landm(i,j,k)==OCEAN).AND.(landm(i+1,j,k)==OCEAN)) THEN
         IF ((landm(i+1,j+1,k)==OCEAN).AND.(landm(i,j+1,k)==OCEAN)) THEN
            CYCLE loop  ! point is ok, surrounded by OCEAN cells, do not change
         END IF
         END IF
         ! if not interupted, then one of the points is not OCEAN, so set co zero
         rowu = ((j-1)*n +(i-1))*2 + 1
         rowv = ((j-1)*n +(i-1))*2 + 2
         colu = rowu*(l) + k
         colv = rowv*(l) + k
         ixBGS%Mzuv%co(colu) = 0
         ixBGS%Mzuv%co(colv) = 0
      END DO loop
      END DO
      END DO
      ixBGS%Mzuv%nnz = nr
      ! make sure Mzuv Mzuv' = I?
!      DO i = 1,nl
!         idx = ixBGS%Mzuv%beg(i)
!         k = SUM(ixBGS%Mzuv%co(idx:idx+l-1))
!         IF (k.EQ.0) WRITE(6,*) "HELP deling door 0", i
!         ixBGS%Mzuv%co(idx:idx+l-1) = ixBGS%Mzuv%co(idx:idx+l-1)/SQRT(REAL(k))
!      ENDDO
      END SUBROUTINE build_Mzuv
    
!**************************************************************
      SUBROUTINE build_Mzuv2(ixBGS)
!     constructs the depth-averaging matrices Mzp and Mzuv
!     such that the matrices are normalized : Mzp*Mzp'=Id, Mzuv*Mzuv'=Id
      USE m_build
      USE m_bgsprec
      USE m_wacsr
!      USE m_wecsr
      USE m_wfree
      USE m_sparsekit
      IMPLICIT none
      INCLUDE 'usr.com'
!     IMPORT/EXPORT
      TYPE (bgsprec),        POINTER :: ixBGS
!     LOCAL
      TYPE (csrmatrix),      POINTER :: ix1, ix2
      INTEGER                        :: nr,nl,ier,i,j,k, dumuv,svs
      INTEGER                        :: rowu, rowv, colu, colv, idx
!     

      nr = ixBGS%Auv%n
      nr = nr   ! real dimension of the problem
      nl = nr/l ! dimension of one layer
      CALL wacsr(nl,nl,ixBGS%Mzuv)
      ixBGS%Mzuv%beg(nl+1) = nl+1
      DO i = 1,nl
         ixBGS%Mzuv%beg(i) = i
         ixBGS%Mzuv%jco(i) = i + (l-1)*nl 
      END DO
      
      ixBGS%Mzuv%co = 1
      ixBGS%Mzuv%nnz = nl
      END SUBROUTINE build_Mzuv2

!**************************************************************
      SUBROUTINE depth_average_block(ixA,ixM,ixN,ixB)
!     constructs the depth-average of the csr-matrix A
!     such that B = MAN
      USE m_build
      USE m_wfree
      USE m_sparsekit
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER :: ixA, ixM, ixN, ixB
!     LOCAL
      TYPE (csrmatrix), POINTER :: ixC  !intermediate matrix
      INTEGER             :: ier
!
      CALL mult(ixA,ixN,ixC,ixA%nnz) ! C = A*N
      CALL mult(ixM,ixC,ixB,ixC%nnz) ! B = M*A*N
      CALL csrfree(ixC)
      END SUBROUTINE depth_average_block
    
!**************************************************************
      SUBROUTINE build_svp2(ixBGS)
!     builds the singular vectors of the pressure
      USE m_bgsprec
      USE m_build
      USE m_csrvec
      USE m_sparsekit
      IMPLICIT none
      INCLUDE 'usr.com'
!     IMPORT/EXPORT
      TYPE (bgsprec),     POINTER :: ixBGS
!     LOCAL
      INTEGER                     :: dim, dim2, i,j,k, count, ier, row
      REAL, DIMENSION(:), POINTER :: testvec,svp5, svp6,scal, ordG
      REAL                        :: factor
      TYPE (csrmatrix),   POINTER :: ixG,ixGT,ixGord
!
      dim = ixBGS%Dw%n
      ALLOCATE(ixBGS%svp1(dim),ixBGS%svp2(dim))
      ixBGS%svp1 = 0
      ixBGS%svp2 = 0
      count = 0
      DO k = 1,l
         DO j = 1,m
            DO i = 1,n
               IF (landm(i,j,k)==OCEAN) THEN
                  count = count+1
                  IF (MOD(i+j,2)==0) THEN
                     ixBGS%svp1(count) = 1
                  ELSE
                     ixBGS%svp2(count) = 1
                  END IF
               END IF
            END DO
         END DO
      END DO
      IF (count/=dim) STOP "dimension of the singular vectors of pressure is wrong"
      
      ! normalize the vectors
      ixBGS%svp1 = ixBGS%svp1/SQRT(sum(ixBGS%svp1))
      ixBGS%svp2 = ixBGS%svp2/SQRT(sum(ixBGS%svp2))
      
      ! test if the singular vectors are indeed singular
      ALLOCATE(testvec(ixBGS%Gw%n))
      testvec = 0 
      CALL csrvec(1D0,ixBGS%Gw,ixBGS%svp1,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp1 is not singular to Gw"
!         STOP "in build_svp svp1 is not singular to Gw"
      END IF
      testvec = 0 
      CALL csrvec(1D0,ixBGS%Gw,ixBGS%svp2,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp2 is not singular to Gw"
!         STOP "in build_svp svp2 is not singular to Gw"
      END IF
      DEALLOCATE(testvec)
      ALLOCATE(testvec(ixBGS%Guv%n))
      testvec = 0
      CALL csrvec(1D0,ixBGS%Guv,ixBGS%svp1,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp1 is not singular to Guv"
!         STOP "in build_svp svp1 is not singular to Guv"
      END IF
      testvec = 0
      CALL csrvec(1D0,ixBGS%Guv,ixBGS%svp2,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "in build_svp svp2 is not singular to Guv"
!         STOP "in build_svp svp2 is not singular to Guv"
      END IF
      DEALLOCATE(testvec)

      ! compute svp3 and svp4 the depth-averaged singular pressure field
      dim = ixBGS%M2Duv%n
      ALLOCATE(ixBGS%svp3(dim),ixBGS%svp4(dim))
      ixBGS%svp3 = 0
      ixBGS%svp4 = 0
      CALL csrvec(1D0,ixBGS%M1,ixBGS%svp1,ixBGS%svp3)
      CALL csrvec(1D0,ixBGS%M1,ixBGS%svp2,ixBGS%svp4)

      ! normalize the vectors
      ixBGS%svp3 = ixBGS%svp3/SQRT(sum(ixBGS%svp3))
      ixBGS%svp4 = ixBGS%svp4/SQRT(sum(ixBGS%svp4))

      ! test if the singular vectors are indeed singular
      ALLOCATE(testvec(ixBGS%GuvM1T%n))
      testvec = 0
      CALL csrvec(1D0,ixBGS%GuvM1T,ixBGS%svp3,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp3 is not singular to GuvMzpT"
         WRITE(*,*) sqrt(dot_product(testvec,testvec))
!         WRITE(*,*) testvec
         !STOP "in build_svp svp3 is not singular to GuvMzpT"
      END IF
      testvec = 0
      CALL csrvec(1D0,ixBGS%GuvM1T,ixBGS%svp4,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp4 is not singular to GuvMzpT"
         WRITE(*,*) sqrt(dot_product(testvec,testvec))
!         WRITE(*,*) testvec
!         STOP "in build_svp svp4 is not singular to GuvMzpT"
      END IF
      DEALLOCATE(testvec)

      END SUBROUTINE build_svp2
      
!**************************************************************
      SUBROUTINE build_svp(ixBGS)
!     builds the singular vectors of the pressure
      USE m_bgsprec
      USE m_build
      USE m_csrvec
      USE m_sparsekit
      USE m_wacsr
      USE m_wfree
      IMPLICIT none
      INCLUDE 'usr.com'
!     IMPORT/EXPORT
      TYPE (bgsprec),     POINTER :: ixBGS
!     LOCAL
      INTEGER                     :: dim, i,j,k, count
      REAL, DIMENSION(:), POINTER :: testvec,scal
      TYPE (csrmatrix),   POINTER :: ixGT
!
      dim = ixBGS%Duv%n
      ALLOCATE(ixBGS%svp1(dim),ixBGS%svp2(dim))
      ALLOCATE(scal(dim))

      ! compute pattern of svp1 and svp2 the singular pressure fields
      dim = ixBGS%Dw%n
      ! ALLOCATE(ixBGS%svp1(dim),ixBGS%svp2(dim))
      ixBGS%svp1 = 0
      ixBGS%svp2 = 0
      count = 0
      DO k = 1,l
         DO j = 1,m
            DO i = 1,n
               IF (landm(i,j,k)==OCEAN) THEN
                  count = count+1
                  IF (MOD(i+j,2)==0) THEN
                     ixBGS%svp1(count) = 1
                  ELSE
                     ixBGS%svp2(count) = 1
                  END IF
               END IF
            END DO
         END DO
      END DO
      IF (count/=dim) STOP "dimension of the singular vectors of pressure is wrong"

      ! normalize the vectors
      ixBGS%svp1 = ixBGS%svp1/SQRT(sum(ixBGS%svp1))
      ixBGS%svp2 = ixBGS%svp2/SQRT(sum(ixBGS%svp2))
      
      ! test if scaling is needed (in case of transposed problem)
      ALLOCATE(testvec(ixBGS%Guv%n))
      testvec = 0 
      CALL csrvec(1D0,ixBGS%Guv,ixBGS%svp1,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         CALL transpose_csr(ixBGS%Guv,ixBGS%Guv%n,ixBGS%Duv%n,ixGT)
         DO i = 1,ixGT%n
            j = ixGT%beg(i)
            k = ixGT%beg(i+1)
            IF (k.GT.j) THEN
               scal(i) = 1D0/abs(ixGT%co(j))
            END IF
         END DO
         CALL csrfree(ixGT)
      ELSE
         scal = 1
      END IF

      DEALLOCATE(testvec)
      DO i = 1,dim
         ixBGS%svp1(i) = scal(i)*ixBGS%svp1(i)
         ixBGS%svp2(i) = scal(i)*ixBGS%svp2(i)
      END DO
      ! normalize the vectors
      ixBGS%svp1 = ixBGS%svp1/SQRT(dot_product(ixBGS%svp1,ixBGS%svp1))
      ixBGS%svp2 = ixBGS%svp2/SQRT(dot_product(ixBGS%svp2,ixBGS%svp2))

      ! test if the singular vectors are indeed singular
      ALLOCATE(testvec(ixBGS%Gw%n))
      testvec = 0 
      CALL csrvec(1D0,ixBGS%Gw,ixBGS%svp1,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
        WRITE(*,*) "WARNING: in build_svp svp1 is not singular to Gw"
      ELSE
!        WRITE(*,*) "in build_svp svp1 is singular to Gw!"
      END IF
      testvec = 0 
      CALL csrvec(1D0,ixBGS%Gw,ixBGS%svp2,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
        WRITE(*,*) "WARNING: in build_svp svp2 is not singular to Gw"
      ELSE
!        WRITE(*,*) "in build_svp svp2 is singular to Gw!"
      END IF
      DEALLOCATE(testvec)
      ALLOCATE(testvec(ixBGS%Guv%n))
      testvec = 0
      CALL csrvec(1D0,ixBGS%Guv,ixBGS%svp1,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
        WRITE(*,*) "WARNING: in build_svp svp1 is not singular to Guv"
      ELSE
!        WRITE(*,*) "in build_svp svp1 is singular to Guv!"
      END IF
      testvec = 0
      CALL csrvec(1D0,ixBGS%Guv,ixBGS%svp2,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
        WRITE(*,*) "WARNING: in build_svp svp2 is not singular to Guv"
      ELSE
!        WRITE(*,*) "in build_svp svp2 is singular to Guv!"
      END IF
      DEALLOCATE(testvec)

      ! compute svp3 and svp4 the depth-averaged singular pressure field
      dim = ixBGS%M2Duv%n
      ALLOCATE(ixBGS%svp3(dim),ixBGS%svp4(dim))
      ixBGS%svp3 = 0
      ixBGS%svp4 = 0
      CALL csrvec(1D0,ixBGS%M1,ixBGS%svp1,ixBGS%svp3)
      CALL csrvec(1D0,ixBGS%M1,ixBGS%svp2,ixBGS%svp4)

      ! normalize the vectors
      ixBGS%svp3 = ixBGS%svp3/SQRT(dot_product(ixBGS%svp3,ixBGS%svp3))
      ixBGS%svp4 = ixBGS%svp4/SQRT(dot_product(ixBGS%svp4,ixBGS%svp4))

      ! test if the singular vectors are indeed singular
      ALLOCATE(testvec(ixBGS%GuvM1T%n))
      testvec = 0
      CALL csrvec(1D0,ixBGS%GuvM1T,ixBGS%svp3,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp3 is not singular to GuvMzpT"
         WRITE(*,*) sqrt(dot_product(testvec,testvec))
         ixBGS%svp3 = 0
      ELSE
!        WRITE(*,*) "in build_svp svp3 is singular to GuvMzpT!"
      END IF
      testvec = 0
      CALL csrvec(1D0,ixBGS%GuvM1T,ixBGS%svp4,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp4 is not singular to GuvMzpT"
         WRITE(*,*) sqrt(dot_product(testvec,testvec))
         ixBGS%svp4 = 0
      ELSE
!        WRITE(*,*) "in build_svp svp4 is singular to GuvMzpT!"
      END IF
      DEALLOCATE(testvec)
      ixBGS%svp1 = 0
      ixBGS%svp2 = 0
      CALL csrvec(1D0,ixBGS%M1T,ixBGS%svp3,ixBGS%svp1)
      CALL csrvec(1D0,ixBGS%M1T,ixBGS%svp4,ixBGS%svp2)

      ! test if the singular vectors are indeed singular
      ALLOCATE(testvec(ixBGS%Gw%n))
      testvec = 0 
      CALL csrvec(1D0,ixBGS%Gw,ixBGS%svp1,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp1 is not singular to Gw"
         ixBGS%svp1 = 0
      ELSE
!        WRITE(*,*) "in build_svp svp1 is singular to Gw!"
      END IF
      testvec = 0 
      CALL csrvec(1D0,ixBGS%Gw,ixBGS%svp2,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp2 is not singular to Gw"
         ixBGS%svp1 = 0
      ELSE
!        WRITE(*,*) "in build_svp svp2 is singular to Gw!"
      END IF
      DEALLOCATE(testvec)
      ALLOCATE(testvec(ixBGS%Guv%n))
      testvec = 0
      CALL csrvec(1D0,ixBGS%Guv,ixBGS%svp1,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp1 is not singular to Guv"
         ixBGS%svp1 = 0
      ELSE
!        WRITE(*,*) "in build_svp svp1 is singular to Guv!"
      END IF
      testvec = 0
      CALL csrvec(1D0,ixBGS%Guv,ixBGS%svp2,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp2 is not singular to Guv"
         ixBGS%svp2 = 0
      ELSE
!        WRITE(*,*) "in build_svp svp2 is singular to Guv!"
      END IF
      DEALLOCATE(testvec,scal)

      END SUBROUTINE build_svp
      
!**************************************************************
      SUBROUTINE build_svp3(ixBGS)
!     builds the singular vectors of the pressure
      USE m_bgsprec
      USE m_build
      USE m_csrvec
      USE m_sparsekit
      USE m_wacsr
      USE m_wfree
      IMPLICIT none
      INCLUDE 'usr.com'
!     IMPORT/EXPORT
      TYPE (bgsprec),     POINTER :: ixBGS
!     LOCAL
      INTEGER                     :: dim, i,j,k, count
      REAL, DIMENSION(:), POINTER :: testvec,testvec2,scal,svp1,svp2
      TYPE (csrmatrix),   POINTER :: ixGT,ixM,ixMT, ixA,ixG
!
      dim = ixBGS%M2Duv%n
      ALLOCATE(ixBGS%svp3(dim),ixBGS%svp4(dim))
      ALLOCATE(scal(dim))
      CALL extract_submatrix(ixBGS%M1T,ixA,ixBGS%M1T%n-dim+1,ixBGS%M1T%n,1,dim)
!      CALL write_sparse_matrix(ixA,'MpzTpart')
      
      ! compute pattern of svp3 and svp4 the singular pressure fields
      ixBGS%svp3 = 0
      ixBGS%svp4 = 0
      count = 0
      DO j = 1,m
         DO i = 1,n
            IF (landm(i,j,l)==OCEAN) THEN
               count = count+1
               IF (MOD(i+j,2)==0) THEN
                  ixBGS%svp3(ixA%jco(count)) = 1
               ELSE
                  ixBGS%svp4(ixA%jco(count)) = 1
               END IF
            END IF
         END DO
      END DO
      IF (count/=dim) STOP "dimension of the singular vectors of pressure is wrong"

      ! normalize the vectors
      ixBGS%svp3 = ixBGS%svp3/SQRT(sum(ixBGS%svp3))
      ixBGS%svp4 = ixBGS%svp4/SQRT(sum(ixBGS%svp4))
      
      ! test if scaling is needed (in case of transposed problem)
      ALLOCATE(testvec(ixBGS%GuvM1T%n),testvec2(ixBGS%GuvM1T%n))
      testvec = 0
      testvec2 = 0
      CALL csrvec(1D0,ixBGS%GuvM1T,ixBGS%svp3,testvec)
      CALL csrvec(1D0,ixBGS%GuvM1T,ixBGS%svp4,testvec2)
      IF ((sqrt(dot_product(testvec,testvec))>1E-12).OR.(sqrt(dot_product(testvec2,testvec2))>1E-12)) THEN
         CALL extract_submatrix(ixBGS%GuvM1T,ixM,2*n*m*(l-1)+1,ixBGS%GuvM1T%n,1,dim)
         CALL transpose_csr(ixM,ixM%n,dim,ixMT) 
!         CALL write_sparse_matrix(ixMT,'MT')
         CALL transpose_csr(ixBGS%GuvM1T,ixBGS%GuvM1T%n,dim,ixGT) 
!         CALL write_sparse_matrix(ixGT,'GT')
         DO i = 1,ixMT%n
            j = ixMT%beg(i)
            k = ixMT%beg(i+1)
            IF (k.GT.j) THEN
               scal(i) = 1D0/abs(ixMT%co(j))
               ixMT%co(j:k-1) = scal(i)*ixMT%co(j:k-1)
               ixGT%co(ixGT%beg(i):ixGT%beg(i+1)-1) = scal(i)*ixGT%co(ixGT%beg(i):ixGT%beg(i+1)-1)
            END IF
         END DO
!         CALL write_sparse_matrix(ixBGS%GuvMzpT,'GuvMzpT')
!         CALL write_sparse_matrix(ixBGS%Guv,'Guv')
!         CALL write_sparse_matrix(ixMT,'MTScal')
!         CALL write_sparse_matrix(ixGT,'GTScal')
!         CALL write_sparse_matrix(ixBGS%MzpT,'MzpT')
         CALL csrfree(ixGT)
         CALL csrfree(ixM)
      ELSE
         scal = 1
      END IF
      DEALLOCATE(testvec)

      DO i = 1,dim
         ixBGS%svp3(i) = scal(i)*ixBGS%svp3(i)
         ixBGS%svp4(i) = scal(i)*ixBGS%svp4(i)
      END DO
      ! normalize the vectors
      ixBGS%svp3= ixBGS%svp3/SQRT(dot_product(ixBGS%svp3,ixBGS%svp3))
      ixBGS%svp4= ixBGS%svp4/SQRT(dot_product(ixBGS%svp4,ixBGS%svp4))

      ! test if the singular vectors are indeed singular
      ALLOCATE(testvec(ixBGS%GuvM1T%n))
      testvec = 0
      CALL csrvec(1D0,ixBGS%GuvM1T,ixBGS%svp3,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp3 is not singular to GuvMzpT"
         WRITE(*,*) sqrt(dot_product(testvec,testvec))
      ELSE
         WRITE(*,*) "in build_svp svp3 is singular to GuvMzpT!"
      END IF
      testvec = 0
      CALL csrvec(1D0,ixBGS%GuvM1T,ixBGS%svp4,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp4 is not singular to GuvMzpT"
         WRITE(*,*) sqrt(dot_product(testvec,testvec))
      ELSE
         WRITE(*,*) "in build_svp svp4 is singular to GuvMzpT!"
      END IF
      DEALLOCATE(testvec)
      dim = ixBGS%Duv%n
      ALLOCATE(ixBGS%svp1(dim),ixBGS%svp2(dim))
      ixBGS%svp1 = 0
      ixBGS%svp2 = 0
      CALL csrvec(1D0,ixBGS%M1T,ixBGS%svp3,ixBGS%svp1)
      CALL csrvec(1D0,ixBGS%M1T,ixBGS%svp4,ixBGS%svp2)

      ! test if the singular vectors are indeed singular
      ALLOCATE(testvec(ixBGS%Gw%n))
      testvec = 0 
      CALL csrvec(1D0,ixBGS%Gw,ixBGS%svp1,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp1 is not singular to Gw"
      ELSE
         WRITE(*,*) "in build_svp svp1 is singular to Gw!"
      END IF
      testvec = 0 
      CALL csrvec(1D0,ixBGS%Gw,ixBGS%svp2,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp2 is not singular to Gw"
      ELSE
         WRITE(*,*) "in build_svp svp2 is singular to Gw!"
      END IF
      DEALLOCATE(testvec)
      ALLOCATE(testvec(ixBGS%Guv%n))
      testvec = 0
      CALL csrvec(1D0,ixBGS%Guv,ixBGS%svp1,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp1 is not singular to Guv"
      ELSE
         WRITE(*,*) "in build_svp svp1 is singular to Guv!"
      END IF
      testvec = 0
      CALL csrvec(1D0,ixBGS%Guv,ixBGS%svp2,testvec)
      IF (sqrt(dot_product(testvec,testvec))>1E-12) THEN
         WRITE(*,*) "WARNING: in build_svp svp2 is not singular to Guv"
      ELSE
         WRITE(*,*) "in build_svp svp2 is singular to Guv!"
      END IF
      DEALLOCATE(testvec,scal)
!      STOP

      END SUBROUTINE build_svp3

!**************************************************************
      SUBROUTINE compute_big_svps(ixBGS,svp1,svp2,dump)
      USE m_build
      USE m_bgsprec
      USE m_sparsekit
      IMPLICIT none
      INCLUDE 'usr.com'
!     IMPORT/EXPORT
      TYPE (bgsprec), POINTER      :: ixBGS
      REAL, DIMENSION(:)           :: svp1,svp2      
      INTEGER                      :: dumw,dump
!     LOCAL
      INTEGER :: nr
      nr = n*m*(l+la)
      svp1=0.0
      svp2=0.0
      svp1(3*nr+dump+1:4*nr) = ixBGS%svp1
      svp2(3*nr+dump+1:4*nr) = ixBGS%svp2
      END SUBROUTINE compute_big_svps       

!**************************************************************
      SUBROUTINE build_diag_blocks(ixBGS)
!     build the diagonal block Aw and the subblocks of Ap
!     Aw is the square part of Dw
!     Ap is the square part of Gw
      USE m_bgsprec
      USE m_sparsekit
      USE m_wfree
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (bgsprec),   POINTER :: ixBGS
!     LOCAL
      TYPE (csrmatrix), POINTER :: ixAwt,ixInvAwt
      INTEGER :: n1,n2
      n1 = ixBGS%Gw%n
      n2 = ixBGS%M1%n
      ! extract Ap
      CALL extract_submatrix(ixBGS%Gw ,ixBGS%Ap ,1,n1,1,n1)
!      CALL extract_submatrix(ixBGS%M1,ixBGS%M11,1,n2,1,n1)
!      CALL extract_submatrix(ixBGS%M1,ixBGS%M12,1,n2,n1+1,n1+n2)
!      CALL transpose_csr(ixBGS%Mzp1,n2,n1,ixBGS%Mzp1T)
!      CALL transpose_csr(ixBGS%Mzp2,n2,n2,ixBGS%Mzp2T)

      ! now extract Aw
      CALL extract_submatrix(ixBGS%Dw,ixBGS%Aw,1,n1,1,n1)

      ! now extract Duv1
      CALL extract_submatrix(ixBGS%Duv,ixBGS%Duv1,1,n1,1,ixBGS%Auv%n)
      ! now extract Guv1
      CALL extract_submatrix(ixBGS%Guv,ixBGS%Guv1,1,ixBGS%Auv%n,1,n1)

      ! build exact inverse of Ap
      CALL build_inverse(ixBGS%Ap,ixBGS%InvAp)

      ! build exact inverse of Aw (via transpose of Aw)
      CALL transpose_csr(ixBGS%Aw,ixBGS%Aw%n,ixBGS%Aw%n,ixAwt)
      CALL build_inverse(ixAwt,ixInvAwt)
      CALL transpose_csr(ixInvAwt,ixInvAwt%n,ixInvAwt%n,ixBGS%InvAw)
      CALL csrfree(ixAwt)
      CALL csrfree(ixInvAwt)
      END SUBROUTINE build_diag_blocks

!**************************************************************
      SUBROUTINE grad_div_stabilization(ixBGS,omega)
!     computes the grad-div stabilized depth-averaged saddle point problem
      USE m_bgsprec
      USE m_sparsekit
      USE m_wfree
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (bgsprec), POINTER       :: ixBGS
      REAL                          :: omega
      
      WRITE(6,*) "Computing grad-div stabilized depth-averaged saddle point matrix"
      ixBGS%omega = omega
      IF (omega == 0) THEN
         CALL copy_csr(ixBGS%MAuv,ixBGS%MAuvGradDiv)
      ELSE
         CALL multadd(ixBGS%MAuv,omega,ixBGS%MGuv,ixBGS%MDuv,ixBGS%MAuvGradDiv)
      END IF
      END SUBROUTINE grad_div_stabilization

!**************************************************************
      SUBROUTINE compute_simpler(ixBGS)
!     computes the blocks needed by the simpler preconditioner
!     for the 'big' saddle point problem
      USE m_bgsprec
      USE m_bgspars
      USE m_sparsekit
      USE m_wfree
      USE m_build
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (bgsprec), POINTER       :: ixBGS
      TYPE (csrmatrix), POINTER     :: ixC
      
      IF (outlev.gt.1) WRITE(6,*) "Extracting 2x2-block-diagonal of Auv"
      CALL extract_block_diagonal(ixBGS%Auv,ixBGS%DAuv,2)
      ! compute exact inverse of 2x2-block diagonal
      IF (outlev.gt.1) WRITE(6,*) "Compute exact inverse of block-diagonal of Auv"
      CALL invert_block_diagonal(ixBGS%DAuv,ixBGS%PDAuv,2)
      ! compute the SIMPLE(R) approximation of the Schur complement
      ! Cp = MpDuv*inv(DMAuv)*GuvMzpT
      IF (outlev.gt.1) WRITE(6,*) "Constructing Modified SIMPLE(R) matrix Cp "
      CALL mult(ixBGS%PDAuv,ixBGS%GuvM1T,ixC,2*ixBGS%GuvM1T%nnz)
      CALL mult(ixBGS%M2Duv,ixC,ixBGS%Cp,2*ixBGS%M2Duv%nnz)
      CALL csrfree(ixC)
      
      END SUBROUTINE compute_simpler

!**************************************************************
      SUBROUTINE compute_precons(ixBGS,type)
!     computes the preconditioners of the diagonal blocks in the
!     BGS preconditioner
      USE m_bgsprec
      USE m_sparsekit
      USE m_bgspars
      ! modules from f90mrilu-libraries
      USE m_cmpprc
      USE m_wacsr
      USE m_waprc
      USE m_iniglb
      USE m_iniprc
      USE m_chkcnt
      USE m_wfree
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (bgsprec),   POINTER :: ixBGS
      INTEGER :: type ! 2 = depth averaged spp, 3 = modified simpler
!     LOCAL
      TYPE (csrmatrix), POINTER :: ixAuv,ixMAuv,ixATS,ixP,ixDMAuv,ixC
!      TYPE (mat), POINTER     :: ixAuv,ixMAuv,ixAS,ixAT, ixP
      INTEGER :: n, blksz, nnzd, nnzoff, nnzlas,ier,i,j
      LOGICAL :: xactelm, gusmod
      REAL    :: epsw, droptol
      REAL    :: time1, time2
      CHARACTER (LEN = 6):: matname

      if (outlev == 3) then
         CALL iniglb(5)
      else
         CALL iniglb(outlev)
      end if

      IF (type.eq.1) THEN
         IF ((prectype=='GD').OR.(prectype=='AC')) THEN
            ! compute mrilu-preconditioner (block-size=1) for MauvGradDiv
            matname = 'MAuvGD'
            CALL copy_csr(ixBGS%MAuvGradDiv,ixMAuv)
         ELSE
            matname = 'MAuv  '
            CALL copy_csr(ixBGS%MAuv,ixMAuv)
         END IF
         IF (outlev.gt.1) WRITE(6,9000) matname
         CALL set_mrilu_pars('a')
         CALL cmpprc(2,ixMAuv,ixBGS%PMAuv)
         CALL chkcnt(ixBGS%PMAuv%mlp,nnzd,nnzoff,nnzlas)
         WRITE(6,9002) matname,int((nnzd + nnzoff + nnzlas)/ixBGS%MAuv%n)
         ! if needed, compute mrilu-preconditioner for MCp
         IF ((prectype=='SI').OR.(prectype=='SR')) THEN
            matname = 'MCp   '
            ! extract 2x2-block diagonal from MAuv
            IF (outlev.gt.1) WRITE(6,*) "Extracting block-diagonal of MAuv"
            CALL extract_block_diagonal(ixBGS%MAuv,ixDMAuv,2)
            ! compute exact inverse of 2x2-block diagonal
            IF (outlev.gt.1) WRITE(6,*) "Compute exact inverse of block-diagonal of MAuv"
            CALL invert_block_diagonal(ixDMAuv,ixBGS%PDMAuv,2)
            ! compute the SIMPLE(R) approximation of the Schur complement
            ! MCp = MDuv*inv(DMAuv)*MGuv
            IF (outlev.gt.1) WRITE(6,*) "Constructing SIMPLE(R) matrix MCp "
            CALL mult(ixBGS%PDMAuv,ixBGS%MGuv,ixC,2*ixBGS%MGuv%nnz)
            CALL mult(ixBGS%MDuv,ixC,ixBGS%MCp,2*ixBGS%MDuv%nnz)
            CALL csrfree(ixC)
            ! compute mrilu-preconditioner (block-size=1) for MCp
            IF (outlev.gt.1) WRITE(6,9000) matname
            CALL copy_csr(ixBGS%MCp,ixC)
            CALL set_mrilu_pars('a')
            CALL cmpprc(1,ixC,ixBGS%PMCp)
            CALL chkcnt(ixBGS%PMCp%mlp,nnzd,nnzoff,nnzlas)
            WRITE(6,9002) matname,int((nnzd + nnzoff + nnzlas)/ixBGS%MCp%n)
         ELSE
            CALL wacsr(1,1,ixBGS%PDMAuv)  ! allocate not used 
            CALL wacsr(1,1,ixBGS%MCp)     ! matrices of bgsprec
            CALL waprc(1,ixBGS%PMCp)      !
         END IF
      ELSEIF (type.eq.2) THEN
         matname = 'Cp'
         ! compute mrilu-preconditioner (block-size=1) for Cp
         IF (outlev.gt.1) WRITE(6,9000) matname
         CALL copy_csr(ixBGS%Cp,ixC)
         CALL set_mrilu_pars('a')
         CALL cmpprc(1,ixC,ixBGS%PCp)
         CALL chkcnt(ixBGS%PCp%mlp,nnzd,nnzoff,nnzlas)
         WRITE(6,9002) matname,int((nnzd + nnzoff + nnzlas)/ixBGS%Cp%n)        
      END IF

      ! compute mrilu-preconditioner (block-size=2) for Auv
      matname = 'Auv   '
      IF (outlev.gt.1) WRITE(6,9000) matname
      CALL copy_csr(ixBGS%Auv,ixAuv)
      CALL set_mrilu_pars('b')
      CALL cmpprc(2,ixAuv,ixBGS%PAuv)
      CALL chkcnt(ixBGS%PAuv%mlp,nnzd,nnzoff,nnzlas)
      WRITE(6,9002) matname, int((nnzd + nnzoff + nnzlas)/ixBGS%Auv%n)

      ! compute mrilu-preconditioner for ATS or the related Schur complement
      IF (xctsc) THEN
         if (type.eq.1) then
            ! extract block-diagonal of Auv
            CALL extract_block_diagonal(ixBGS%Auv,ixBGS%DAuv,2)
            CALL invert_block_diagonal(ixBGS%DAuv,ixBGS%PDAuv,2)
         end if
         matname = 'ATS+..'
         CALL cpu_time(time1)
         CALL cmp_schurcomp_ATS(ixBGS,ixATS)
         CALL cpu_time(time2)
         IF (outlev.gt.1) WRITE(6,*) 'Time to construct the ATS-Schur complement: ',time2-time1
      ELSE
         matname = 'ATS   '
         CALL copy_csr(ixBGS%ATS,ixATS)
      END IF
      ! START positive diagonal test
      j = 0
!      CALL extract_block_diagonal(ixBGS%ATS,ixC,1)
      ! check diagonal
!      IF ((ixC%n.eq.ixBGS%ATS%n).AND.(ixC%nnz.eq.ixC%n)) THEN
!         DO i = 1,ixC%n
!            IF (ixC%co(i).LT.0D0) THEN
!               WRITE(6,*) "found negative diagonal entry in ATS", i,ixC%jco(i)
!               j = j+1
!            END IF
!         END DO
!      ELSE
!         WRITE (6,*) "error in diagonal of ATS" 
!      END IF
!      WRITE(6,*) "found total of ",j, "negative diagonal entries in ATS " 
      ! END positive diagonal test
!      CALL csrfree(ixC)
      IF (outlev.gt.1) WRITE(6,9000) matname
      CALL set_mrilu_pars('c')
      CALL cmpprc(1,ixATS,ixBGS%PATS) ! block size = 1
      CALL chkcnt(ixBGS%PATS%mlp,nnzd,nnzoff,nnzlas)
      WRITE(6,9002) matname, int((nnzd + nnzoff + nnzlas)/ixBGS%ATS%n)

9000  FORMAT (" Constructing preconditioner for ",a6)
9002  FORMAT (" nnz/row in the MRILU-decompostion of ",a6,":",I4)   

      END SUBROUTINE compute_precons

!**************************************************************
      SUBROUTINE cmp_schurcomp_ATS(ixBGS,ixATS)
        ! computes a very good approximation of the Schur complement ATS
        USE m_build
        USE m_bgsprec
        USE m_sparsekit
        USE m_wfree
        IMPLICIT none
        ! IMPORT/EXPORT
        TYPE (bgsprec),   POINTER :: ixBGS
        TYPE (csrmatrix), POINTER :: ixATS
        ! LOCAL
        TYPE (csrmatrix), POINTER :: ixA,ixB,ixC,ixD
        INTEGER              :: nnzA
        !
        IF ((ixBGS%BwTS%nnz==0).OR.((ixBGS%BTSuv%nnz==0).AND.(ixBGS%BTSw%nnz==0))) THEN
           CALL copy_csr(ixBGS%ATS,ixATS)
           nnzA = ixBGS%ATS%nnz
        ELSE
           ! first compute Ap\BwTS
           CALL mult(ixBGS%InvAp,ixBGS%BwTS,ixA,4*ixBGS%InvAp%nnz) ! A = Gw1\BwTS
           ! compute Guv*(Ap\BwTS)
           CALL mult(ixBGS%Guv1,ixA,ixB,10*ixA%nnz)    ! B = Guv*A
           ! compute approximation of Auv\(Guv1*(Ap\BwTS))
           CALL mult(ixBGS%PDAuv,ixB,ixC,4*ixB%nnz)   ! C = DAuv\B
           CALL csrfree(ixA); CALL csrfree(ixB)       ! delete A and B
           ! compute Aw\(Duv1*C)
           CALL mult(ixBGS%Duv1,ixC,ixA,2*ixC%nnz)
           CALL mult(ixBGS%InvAw,ixA,ixB,4*ixA%nnz)
           ! now compute ATS + BTSuv*C - BTSw*B
           CALL multadd(ixBGS%ATS,1D0,ixBGS%BTSuv,ixC,ixD)
           CALL multadd(ixD,-1D0,ixBGS%BTSw,ixB,ixATS)
           CALL csrfree(ixA); CALL csrfree(ixB)
           CALL csrfree(ixC); CALL csrfree(ixD) 
           nnzA =        ixBGS%BwTS%nnz  + ixBGS%Ap%nnz
           nnzA = nnzA + ixBGS%Guv1%nnz
           nnzA = nnzA + ixBGS%PDAuv%nnz + ixBGS%Duv1%nnz
           nnzA = nnzA + ixBGS%Aw%nnz    + ixBGS%ATS%nnz
           nnzA = nnzA + ixBGS%BTSuv%nnz + ixBGS%BTSw%nnz
        END IF
        WRITE(6,*) "nnz/row ATS      :", ixBGS%ATS%nnz/ixATS%n 
        WRITE(6,*) "nnz/row ATS2     :", ixATS%nnz/ixATS%n
        WRITE(6,*) "nnz/row applATS2 :", nnzA/ixATS%n
!       CALL write_sparse_matrix(ixATS,'STS')
      END SUBROUTINE cmp_schurcomp_ATS

!**************************************************************
      SUBROUTINE set_mrilu_pars(system)
!     sets the mrilu parameters
      USE m_bgspars
      USE m_mrilupars
      USE m_iniprc
      IMPLICIT none
!     IMPORT/EXPORT
      CHARACTER   :: system
!
      LOGICAL     :: a,b,c,d,t,u
      REAL        :: e,f,g,h,i,j,k,l,m,n,p,q,r,s
      INTEGER     :: o 
!
      IF (system == 'a') THEN
         a = MAuv_cutmck;  b = MAuv_scarow;   c = MAuv_xactelm
         d = MAuv_clsonce; e = MAuv_nlsfctr;  f = MAuv_epsw
         g = MAuv_elmfctr; h = MAuv_gusfctr;  i = MAuv_redfctr;
         j = MAuv_schtol;  k = MAuv_denslim;  l = MAuv_globfrac;
         m = MAuv_locfrac; n = MAuv_sparslim; o = MAuv_ilutype;
         p = MAuv_droptol; q = MAuv_compfct;  r = MAuv_cpivtol;
         s = MAuv_lutol;   t = MAuv_singlu;   u = .false.
      ELSEIF (system == 'b') THEN
         a = Auv_cutmck;  b = Auv_scarow;   c = Auv_xactelm
         d = Auv_clsonce; e = Auv_nlsfctr;  f = Auv_epsw
         g = Auv_elmfctr; h = Auv_gusfctr;  i = Auv_redfctr;
         j = Auv_schtol;  k = Auv_denslim;  l = Auv_globfrac;
         m = Auv_locfrac; n = Auv_sparslim; o = Auv_ilutype;
         p = Auv_droptol; q = Auv_compfct;  r = Auv_cpivtol;
         s = Auv_lutol;   t = Auv_singlu;   u = .false.
      ELSEIF (system == 'c') THEN
         a = ATS_cutmck;  b = ATS_scarow;   c = ATS_xactelm
         d = ATS_clsonce; e = ATS_nlsfctr;  f = ATS_epsw
         g = ATS_elmfctr; h = ATS_gusfctr;  i = ATS_redfctr;
         j = ATS_schtol;  k = ATS_denslim;  l = ATS_globfrac;
         m = ATS_locfrac; n = ATS_sparslim; o = ATS_ilutype;
         p = ATS_droptol; q = ATS_compfct;  r = ATS_cpivtol;
         s = ATS_lutol;   t = ATS_singlu;   u = .false.
      ELSEIF (system == 'd') THEN
         a = cutmck;  b = scarow;   c = xactelm
         d = clsonce; e = nlsfctr;  f = epsw
         g = elmfctr; h = gusfctr;  i = redfctr;
         j = schtol;  k = denslim;  l = globfrac;
         m = locfrac; n = sparslim; o = ilutype;
         p = droptol; q = compfctr; r = cpivtol;
         s = lutol;   t = singlu;   u = .false.
      ELSE
         WRITE(6,*) "mrilupars not specified, using standard values"
      END IF
      CALL iniprc (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) 
      END SUBROUTINE set_mrilu_pars

!**************************************************************
      SUBROUTINE apply_bilu(ixBGS,b,ord,variant,spptype)
!     computes the action of the inverse of the BGS preconditioner
!     on the vector b, in fact the equation BGS*x = b is solved
      ! modules from bgs library
      USE m_bgsprec
      !
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (bgsprec), POINTER                   :: ixBGS
      REAL, DIMENSION(:), INTENT(INOUT), TARGET :: b
      INTEGER, DIMENSION(:), INTENT(IN)         :: ord
      INTEGER :: variant ! schurcomplement: 1 = ST, 2 = uv, 3 = p/w 
      INTEGER :: spptype ! 1 = depth averaged spp, 2 = modified simpler
!     LOCAL
!      INTEGER :: var

      IF (variant.lt.4) THEN
         var = variant
      ELSE
         var = MOD(var+1,3)+1
!         WRITE(6,*) "applying variant:", var
      END IF
      
      IF (var == 1) THEN
         CALL apply_bilu1(ixBGS,b,ord,spptype)
      ELSEIF (var == 2) THEN
         CALL apply_bilu2(ixBGS,b,ord,spptype)
      ELSEIF ((var == 3).AND.(.NOT.ixBGS%transposed)) THEN
         CALL apply_bilu3(ixBGS,b,ord,spptype)
      ELSEIF ((var == 3).AND.ixBGS%transposed) THEN
         CALL apply_bilu3 (ixBGS,b,ord,spptype)
      ELSE
         WRITE(6,*) "Warning variant in apply_bilu (bgskit) not defined"
         continue
      END IF

      END SUBROUTINE apply_bilu

!**************************************************************
      SUBROUTINE apply_bilu1(ixBGS,b,ord,type)
!     computes the action of the inverse of the BGS preconditioner
!     on the vector b, in fact the equation BGS*x = b is solved
      ! modules from mrilu-libraries
      USE m_build
      USE m_csrvec
      USE m_inisol
      USE m_iniglb
      USE m_applprc
      USE m_solprc
      ! modules from bgs library
      USE m_bgsprec
      USE m_bgspars
      USE m_sparsekit
      USE m_bgsgmres
      !
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (bgsprec), POINTER                   :: ixBGS
      REAL, DIMENSION(:), INTENT(INOUT), TARGET :: b
      INTEGER, DIMENSION(:), INTENT(IN)         :: ord
      INTEGER :: type ! 1 = depth averaged spp, 2 = modified simpler
!     LOCAL
      TYPE (prcmatrix),   POINTER :: Prc
      REAL, DIMENSION(:), POINTER :: x
      REAL, DIMENSION(:), POINTER :: xuv,xw,xp,xTS,xdumw,xdump,xz
      REAL, DIMENSION(:), POINTER :: buv,bw,bp,bTS,bdumw,bdump
      REAL, DIMENSION(:), POINTER :: yuv,yuv2,yw,yw2,yp,yTS,yTS2,yzp,yzuv,xzuv,yz
      REAL, DIMENSION(:), POINTER :: buvzp, yuvzp ! dummy arrays
      REAL, DIMENSION(:), POINTER :: xzp, xtilp, bzp, bzuv
      INTEGER                     :: nuv,nw,np,nTS,ndumw,ndump,nzp,nzuv, nz
      INTEGER                     :: n,i
      INTEGER, DIMENSION(:), ALLOCATABLE :: ord2
      REAL :: nrm
  
      call iniglb(outlev)
      ! first reorder rhs 
      b = b(ord)
      ! extract dimensions
      nuv   = ixBGS%PAuv%n;     nw    = ixBGS%Gw%n;      np  = ixBGS%Dw%n;
      nTS   = ixBGS%PATS%n;
      ndumw = ixBGS%Adumw%n;    ndump = ixBGS%Adump%n;  
      nzp   = ixBGS%M1%n;      
      if (type.eq.1) then
         nzuv  = ixBGS%Mzuv%n 
         nz    = nzuv+nzp
      else
         nzuv  = 0
         nz    = nuv+nzp
      end if
      ! set pointers
      n = size(b,1)
      ALLOCATE(x(n))
      xuv   => x(1:nuv);       buv   => b(1:nuv);       n = nuv
      xdumw => x(n+1:n+ndumw); bdumw => b(n+1:n+ndumw); n = n+ndumw
      xw    => x(n+1:n+nw);    bw    => b(n+1:n+nw);    n = n+nw
      xdump => x(n+1:n+ndump); bdump => b(n+1:n+ndump); n = n+ndump
      xp    => x(n+1:n+np);    bp    => b(n+1:n+np);    n = n+np
      xTS   => x(n+1:n+nTS);   bTS   => b(n+1:n+nTS);   n = n+nTS
      
      bp = -bp
      x = 0
      ! allocate the dummy arrays
      ALLOCATE(yp(np),yw(nw),yw2(nw),yzp(nzp),bzp(nzp),bzuv(nuv))
      ALLOCATE(yzuv(nuv),yuv(nuv),yTS(nTS),yTS2(nTS),yz(nz),yuv2(nuv))
      ALLOCATE(xtilp(np),xzp(nzp),xzuv(nuv),xz(nz))
      ALLOCATE(yuvzp(nz),buvzp(nz))

      ! a) solve dummy equations
      CALL dsolve(ixBGS%Adump,bdump,xdump) 
      CALL dsolve(ixBGS%Adumw,bdumw,xdumw)
      yp = 0
      CALL csrvec(1D0,ixBGS%Bpdumw,xdumw,yp)
      bp = bp-yp

      ! b) solve the transformed pressure xp
      CALL usolve(ixBGS%Ap,bw,xp(1:nw))

      ! c) and d) solve saddle point problem
      bzp = 0
      CALL csrvec(1D0,ixBGS%M2,bp,bzp)          ! bzp = M2*bp
      yuv = buv
      CALL csrvec(-1D0,ixBGS%Guv1,xp(1:nw),yuv) ! yuv = buv-Guv*xtilp
      IF (type.eq.1) THEN
         ! c) solve depth-averaged saddle point problem in order
         !    to obtain the depth-averaged pressure
         !    we use the Artficial Compressibillity preconditioner
         !    and are only interested in the pressure part of the solution
         
         !    first construct the right-hand-side of the depth-averaged equations
         bzuv = 0
         CALL csrvec(1D0,ixBGS%Mzuv,yuv,bzuv)  ! bzuv = Mzuv*(buv-Guv*xtilp)
         !    then solve equations via gmres with spp preconditioner
         CALL inisol(cgtype,mgmres,maxnits,spploctolabs,spploctolred)
         yz = 0
         xz(1:nzuv) = bzuv
         xz(nzuv+1:nz) = bzp
         CALL sppgmres(mgmres,nz,ixBGS,yz,xz,prectype)
         xzp = yz(nzuv+1:nz)
         
         ! correct depth-averaged pressure for constant field
!         xzp = xzp - dot_product(xzp,ixBGS%svp3)*ixBGS%svp3
!         xzp = xzp - dot_product(xzp,ixBGS%svp4)*ixBGS%svp4
         
         ! d) solve the velocity field
         CALL csrvec(-1D0,ixBGS%GuvM1T,xzp,yuv)
         ! yuv = buv - Guv*xtilp - GuvM1T*xzp
         CALL inisol(cgtype,mgmres,maxnits,loctolabs,loctolred)
         Prc => ixBGS%PAuv
         xuv = 0
         IF (locgmres) THEN ! solve Auv iteratively with preconditioned gmres
            CALL prcgmres(mgmres,nuv,ixBGS%Auv,Prc,Prc%mlp%perm,xuv,yuv)
         ELSE               ! apply preconditioner of Auv once
            CALL applprc(Prc,xuv,yuv)
         END IF
      ELSEIF (type.eq.2) THEN
         ! c) and d) solve big Saddle Point Problem of
         !    depth averaged pressure and velocity field using a
         !    Simpler preconditioner
         yuvzp = 0
         buvzp(1:nuv) = yuv
         buvzp(nuv+1:nuv+nzp) = bzp
         IF (locgmres) THEN
!            CALL inisol(cgtype,mgmres,maxnits,spploctolabs,spploctolred)
            CALL inisol(cgtype,mgmres,maxnits,1D-6,1D-6)
            CALL bigsppfgmres(mgmres,nuv+nzp,ixBGS,yuvzp,buvzp)
         ELSE
            yuvzp = buvzp
            CALL apply_bigsppprec (ixbgs, yuvzp, 'SR')
         END IF
         xuv = yuvzp(1:nuv)
         xzp = yuvzp(nuv+1:nuv+nzp)
         ! correct pressure for constant field
!         xzp = xzp - dot_product(xzp,ixBGS%svp3)*ixBGS%svp3
!         xzp = xzp - dot_product(xzp,ixBGS%svp4)*ixBGS%svp4         
      END IF
      xp(nw+1:np) = xzp
      
      ! e) solve the vertical velocity-field
      yw = bp(1:nw)
      CALL csrvec(-1D0,ixBGS%Duv1,xuv,yw)
      CALL lsolve(ixBGS%Aw,yw,xw)
      
      ! f and g) solve the temperature and salinity equation
      yTS = bTS
      CALL csrvec(-1D0,ixBGS%BTSuv,xuv,yTS)
      CALL csrvec(-1D0,ixBGS%BTSw ,xw ,yTS)
      Prc => ixBGS%PATS
      CALL inisol(cgtype,mgmres,mgmres,loctolabs,loctolred)
      IF (locgmres) THEN ! solve ATS iteratively with preconditioned gmres
         IF (ilu) THEN
            CALL SCfgmres(mgmres,nTS,ixBGS,xTS,yTS,1,type)
         ELSE
            CALL prcfgmres(mgmres,nTS,ixBGS%ATS,Prc,Prc%mlp%perm,xTS,yTS)
         END IF
      ELSE               ! apply preconditioner of ATS once
         CALL applprc(Prc,xTS,yTS)
      END IF

      IF (ilu) THEN
         ! compute correction (in fact application of inverse of U)
         yw = 0; yzp = 0; yp = 0; yuv = 0; yuv2 = 0;
         ! h) compute correction for xp
         CALL csrvec(1D0,ixBGS%BwTS,xTS,yw)         ! yw = BwTS*xTS
         CALL usolve(ixBGS%Ap,yw,yp(1:nw))          ! yp = Ap\(BwTS*xTS)
         ! i) compute correction for xuv and xzp
         CALL csrvec(-1D0,ixBGS%Guv1,yp(1:nw),yuv2)
         CALL csrvec(1D0,ixBGS%BuvTS,xTS,yuv2)      ! yuv2 = BuvTS*xTS-Guv1*yp
         IF (type.eq.1) THEN
            CALL applprc(ixBGS%PAuv,yuv,yuv2)       ! yuv = PAuv\yuv2
            yzp = 0
         ELSEIF (type.eq.2) THEN
            buvzp = 0; buvzp(1:nuv) = yuv2;
            !         CALL inisol(cgtype,mgmres,maxnits,spploctolabs,spploctolred)
            !         CALL inisol(cgtype,mgmres,maxnits,1D-6,1D-6)
            !         CALL bigsppfgmres(mgmres,nuv+nzp,ixBGS,yuvzp,buvzp) ! solve bigSPP
            yuvzp = buvzp
            CALL apply_bigsppprec (ixbgs, yuvzp, 'SR')
            yuv = yuvzp(1:nuv)                      ! yuv = SPP\[yuv2;0]
            yzp = yuvzp(nuv+1:nuv+nzp)
            ! correct pressure for constant field
!            yzp = yzp - dot_product(yzp,ixBGS%svp3)*ixBGS%svp3
!            yzp = yzp - dot_product(yzp,ixBGS%svp4)*ixBGS%svp4
         END IF
         ! j) compute correction for xw
         yw = 0; yw2 = 0;
         CALL csrvec(-1D0,ixBGS%Duv1,yuv,yw2)
         CALL lsolve(ixBGS%Aw,yw2,yw);              ! yw = -Aw\(Duv1*yuv)
         ! k) add corrections to x
         xp(1:nw) = xp(1:nw) - yp(1:nw);
         xuv = xuv - yuv;
         xp(nw+1:np) = xp(nw+1:np) - yzp;
         xw  = xw  - yw;
      END IF
      ! l) compute transformation of the pressure
      xzp = xp(nw+1:np); xp(nw+1:np) = 0
      CALL csrvec(1D0,ixBGS%M1T,xzp,xp);
      xp = xp - dot_product(xp,ixBGS%svp1)*ixBGS%svp1
      xp = xp - dot_product(xp,ixBGS%svp2)*ixBGS%svp2         

      ! finally reorder solution and put it in b
      b(ord) = x

      DEALLOCATE(x)
      DEALLOCATE(yp,yw,yw2,yzp,bzp,bzuv)
      DEALLOCATE(yzuv,yuv,yuv2,yTS,yTS2,yz)
      DEALLOCATE(xtilp,xzp,xzuv,xz)
      DEALLOCATE(yuvzp,buvzp)
      CALL inisol(cgtype,mgmres,maxnits,tolred,tolabs)
 !     CALL iniglb(2)
      END SUBROUTINE apply_bilu1

!**************************************************************
      SUBROUTINE apply_bilu2(ixBGS,b,ord,type)
!     computes the action of the inverse of the BGS preconditioner
!     on the vector b, in fact the equation BGS*x = b is solved
      ! modules from mrilu-libraries
      USE m_build
      USE m_csrvec
      USE m_inisol
      USE m_iniglb
      USE m_applprc
      USE m_solprc
      ! modules from bgs library
      USE m_bgsprec
      USE m_bgspars
      USE m_sparsekit
      USE m_bgsgmres
      !
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (bgsprec), POINTER                   :: ixBGS
      REAL, DIMENSION(:), INTENT(INOUT), TARGET :: b
      INTEGER, DIMENSION(:), INTENT(IN)         :: ord
      INTEGER :: type ! 1 = depth averaged spp, 2 = modified simpler
!     LOCAL
      TYPE (prcmatrix),   POINTER :: Prc
      REAL, DIMENSION(:), POINTER :: x
      REAL, DIMENSION(:), POINTER :: xuv,xw,xp,xTS,xdumw,xdump,xz
      REAL, DIMENSION(:), POINTER :: buv,bw,bp,bTS,bdumw,bdump
      REAL, DIMENSION(:), POINTER :: yuv,yuv2,yw,yw2,yp,yTS,yTS2,yzp,yzuv,xzuv,yz
      REAL, DIMENSION(:), POINTER :: buvzp, yuvzp ! dummy arrays
      REAL, DIMENSION(:), POINTER :: xzp, xtilp, bzp, bzuv
      INTEGER                     :: nuv,nw,np,nTS,ndumw,ndump,nzp,nzuv, nz
      INTEGER                     :: n,i
      INTEGER, DIMENSION(:), ALLOCATABLE :: ord2
      REAL :: nrm
  
      call iniglb(outlev)
      ! first reorder rhs 
      b = b(ord)
      ! extract dimensions
      nuv   = ixBGS%PAuv%n;     nw    = ixBGS%Gw%n;      np  = ixBGS%Dw%n;
      nTS   = ixBGS%PATS%n;
      ndumw = ixBGS%Adumw%n;    ndump = ixBGS%Adump%n;  
      nzp   = ixBGS%M1%n;      
      if (type.eq.1) then
         nzuv  = ixBGS%Mzuv%n 
         nz    = nzuv+nzp
      else
         nzuv  = 0
         nz    = nuv+nzp
      end if
      ! set pointers
      n = size(b,1)
      ALLOCATE(x(n))
      xuv   => x(1:nuv);       buv   => b(1:nuv);       n = nuv
      xdumw => x(n+1:n+ndumw); bdumw => b(n+1:n+ndumw); n = n+ndumw
      xw    => x(n+1:n+nw);    bw    => b(n+1:n+nw);    n = n+nw
      xdump => x(n+1:n+ndump); bdump => b(n+1:n+ndump); n = n+ndump
      xp    => x(n+1:n+np);    bp    => b(n+1:n+np);    n = n+np
      xTS   => x(n+1:n+nTS);   bTS   => b(n+1:n+nTS);   n = n+nTS
      
      bp = -bp
      x = 0
      ! allocate the dummy arrays
      ALLOCATE(yp(np),yw(nw),yw2(nw),yzp(nzp),bzp(nzp),bzuv(nuv))
      ALLOCATE(yzuv(nuv),yuv(nuv),yTS(nTS),yTS2(nTS),yz(nz),yuv2(nuv))
      ALLOCATE(xtilp(np),xzp(nzp),xzuv(nuv),xz(nz))
      ALLOCATE(yuvzp(nz),buvzp(nz))

      ! a) solve dummy equations
      CALL dsolve(ixBGS%Adump,bdump,xdump) 
      CALL dsolve(ixBGS%Adumw,bdumw,xdumw)
      yp = 0
      CALL csrvec(1D0,ixBGS%Bpdumw,xdumw,yp)
      bp = bp-yp
      
      ! b) solve the vertical velocity-field
      CALL lsolve(ixBGS%Aw,bp(1:nw),xw)
      
      ! c and d) solve the temperature and salinity equation
      yTS = bTS
      CALL csrvec(-1D0,ixBGS%BTSw ,xw ,yTS)
!      CALL applprc(ixBGS%PATS,xTS,yTS)
      IF (locgmres) THEN ! solve ATS iteratively with preconditioned gmres
         Prc => ixBGS%PATS
         CALL inisol(cgtype,mgmres,mgmres,loctolabs,loctolred)
         CALL prcfgmres(mgmres,nTS,ixBGS%ATS,Prc,Prc%mlp%perm,xTS,yTS)
      ELSE               ! apply preconditioner of ATS once
         CALL applprc(ixBGS%PATS,xTS,yTS)
      END IF

      ! e) solve the transformed pressure xp
      yw = bw
      xp(1:nw) = 0
      CALL csrvec(-1D0,ixBGS%BwTS,xTS,yw)
      CALL usolve(ixBGS%Ap,yw,xp(1:nw))

      ! f) and g) solve saddle point problem
      bzp = 0
      CALL csrvec(1D0,ixBGS%M2,bp,bzp)       ! bzp = M2*bp
      yuv = buv
      CALL csrvec(-1D0,ixBGS%BuvTS,xTS,yuv)    
      CALL csrvec(-1D0,ixBGS%Guv1,xp(1:nw),yuv)
      buvzp(1:nuv) = yuv
      buvzp(nuv+1:nuv+nzp) = bzp
      yuvzp = 0
      IF (locgmres) THEN ! solve SC iteratively with preconditioned gmres
         IF (ilu) THEN
            CALL SCgmres(mgmres,nz,ixBGS,yuvzp,buvzp,2,type)
         ELSE
            CALL bigsppfgmres(mgmres,nz,ixBGS,yuvzp,buvzp)
         END IF
      ELSE               ! apply preconditioner of ATS once
         yuvzp = buvzp
         CALL apply_bigsppprec (ixbgs, yuvzp , 'SR')
      END IF

      xuv = yuvzp(1:nuv)
      xzp = yuvzp(nuv+1:nuv+nzp)
      ! correct pressure for constant field
!      xzp = xzp - dot_product(xzp,ixBGS%svp3)*ixBGS%svp3
!      xzp = xzp - dot_product(xzp,ixBGS%svp4)*ixBGS%svp4         
      xp(nw+1:np) = xzp

      IF (ilu) THEN
         ! compute correction (in fact application of inverse of U)
         yw = 0; yzp = 0; yTS = 0; yTS2 = 0; yp = 0; yw2 = 0;
         ! h) compute correction for xw
         CALL csrvec(1D0,ixBGS%Duv1,xuv,yw2)
         CALL lsolve(ixBGS%Aw,yw2,yw);              ! yw = Aw\(Duv1*xuv)
         ! i) compute correction for xTS
         CALL csrvec( 1D0,ixBGS%BTSuv,xuv,yTS2)
         CALL csrvec(-1D0,ixBGS%BTSw, yw, yTS2)
         CALL applprc(ixBGS%PATS,yTS,yTS2)          ! yTS = ATS\(BTSuv*xuv - BTSw*yw)
         ! j) compute correction for xp
         yw2 = 0;
         CALL csrvec(1D0,ixBGS%BwTS,yTS,yw2);
         CALL usolve(ixBGS%Ap,yw2,yp(1:nw));        ! yp = Ap\(BwTS*yTS)
         ! k) add corrections to x
         xw = xw - yw;
         xTS = xTS - yTS;
         xp(1:nw) = xp(1:nw) + yp(1:nw);
      END IF
      ! l) compute transformation of the pressure
      xzp = xp(nw+1:np); xp(nw+1:np) = 0
      CALL csrvec(1D0,ixBGS%M1T,xzp,xp);
      xp = xp - dot_product(xp,ixBGS%svp1)*ixBGS%svp1
      xp = xp - dot_product(xp,ixBGS%svp2)*ixBGS%svp2         

      ! finally reorder solution and put it in b
      b(ord) = x

      DEALLOCATE(x)
      DEALLOCATE(yp,yw,yw2,yzp,bzp,bzuv)
      DEALLOCATE(yzuv,yuv,yuv2,yTS,yTS2,yz)
      DEALLOCATE(xtilp,xzp,xzuv,xz)
      DEALLOCATE(yuvzp,buvzp)
      CALL inisol(cgtype,mgmres,maxnits,tolred,tolabs)
 !     CALL iniglb(2)
      END SUBROUTINE apply_bilu2

!**************************************************************
      SUBROUTINE apply_bilu3(ixBGS,b,ord,type)
!     computes the action of the inverse of the BGS preconditioner
!     on the vector b, in fact the equation BGS*x = b is solved
      ! modules from mrilu-libraries
      USE m_build
      USE m_csrvec
      USE m_inisol
      USE m_iniglb
      USE m_applprc
      USE m_solprc
      ! modules from bgs library
      USE m_bgsprec
      USE m_bgspars
      USE m_sparsekit
      USE m_bgsgmres
      !
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (bgsprec), POINTER                   :: ixBGS
      REAL, DIMENSION(:), INTENT(INOUT), TARGET :: b
      INTEGER, DIMENSION(:), INTENT(IN)         :: ord
      INTEGER :: type ! 1 = depth averaged spp, 2 = modified simpler
!     LOCAL
      TYPE (prcmatrix),   POINTER :: Prc
      REAL, DIMENSION(:), POINTER :: x
      REAL, DIMENSION(:), POINTER :: xuv,xw,xp,xTS,xdumw,xdump,xz
      REAL, DIMENSION(:), POINTER :: buv,bw,bp,bTS,bdumw,bdump
      REAL, DIMENSION(:), POINTER :: yuv,yuv2,yw,yw2,yp,yTS,yTS2,yzp,yzuv,xzuv,yz
      REAL, DIMENSION(:), POINTER :: buvzp, yuvzp ! dummy arrays
      REAL, DIMENSION(:), POINTER :: xzp, xtilp, bzp, bzuv
      INTEGER                     :: nuv,nw,np,nTS,ndumw,ndump,nzp,nzuv, nz
      INTEGER                     :: n,i
      INTEGER, DIMENSION(:), ALLOCATABLE :: ord2
      REAL :: nrm
  
      call iniglb(outlev)
      ! first reorder rhs 
      b = b(ord)
      ! extract dimensions
      nuv   = ixBGS%PAuv%n;     nw    = ixBGS%Gw%n;      np  = ixBGS%Dw%n;
      nTS   = ixBGS%PATS%n;
      ndumw = ixBGS%Adumw%n;    ndump = ixBGS%Adump%n;  
      nzp   = ixBGS%M1%n;      
      if (type.eq.1) then
         nzuv  = ixBGS%Mzuv%n 
         nz    = nzuv+nzp
      else
         nzuv  = 0
         nz    = nuv+nzp
      end if
      ! set pointers
      n = size(b,1)
      ALLOCATE(x(n))
      xuv   => x(1:nuv);       buv   => b(1:nuv);       n = nuv
      xdumw => x(n+1:n+ndumw); bdumw => b(n+1:n+ndumw); n = n+ndumw
      xw    => x(n+1:n+nw);    bw    => b(n+1:n+nw);    n = n+nw
      xdump => x(n+1:n+ndump); bdump => b(n+1:n+ndump); n = n+ndump
      xp    => x(n+1:n+np);    bp    => b(n+1:n+np);    n = n+np
      xTS   => x(n+1:n+nTS);   bTS   => b(n+1:n+nTS);   n = n+nTS
      
      bp = -bp
      x = 0
      ! allocate the dummy arrays
      ALLOCATE(yp(np),yw(nw),yw2(nw),yzp(nzp),bzp(nzp),bzuv(nuv))
      ALLOCATE(yzuv(nuv),yuv(nuv),yTS(nTS),yTS2(nTS),yz(nz),yuv2(nuv))
      ALLOCATE(xtilp(np),xzp(nzp),xzuv(nuv),xz(nz))
      ALLOCATE(yuvzp(nz),buvzp(nz))

      ! a) solve dummy equations
      CALL dsolve(ixBGS%Adump,bdump,xdump) 
      CALL dsolve(ixBGS%Adumw,bdumw,xdumw)
      yp = 0
      CALL csrvec(1D0,ixBGS%Bpdumw,xdumw,yp)
      bp = bp-yp

      ! b) and c) solve saddle point problem
      bzp = 0
      CALL csrvec(1D0,ixBGS%M2,bp,bzp)          ! bzp = M2*bp
      yuv = buv
      IF (type.eq.1) THEN
         ! c) solve depth-averaged saddle point problem in order
         !    to obtain the depth-averaged pressure
         !    we use the Artficial Compressibillity preconditioner
         !    and are only interested in the pressure part of the solution
         
         !    first construct the right-hand-side of the depth-averaged equations
         bzuv = 0
         CALL csrvec(1D0,ixBGS%Mzuv,yuv,bzuv)  ! bzuv = Mzuv*(buv-Guv*xtilp)
         !    then solve equations via gmres with spp preconditioner
         CALL inisol(cgtype,mgmres,maxnits,spploctolabs,spploctolred)
         yz = 0
         xz(1:nzuv) = bzuv
         xz(nzuv+1:nz) = bzp
         CALL sppgmres(mgmres,nz,ixBGS,yz,xz,prectype)
         xzp = yz(nzuv+1:nz)
         
         ! correct depth-averaged pressure for constant field
!         xzp = xzp - dot_product(xzp,ixBGS%svp3)*ixBGS%svp3
!         xzp = xzp - dot_product(xzp,ixBGS%svp4)*ixBGS%svp4
         
         ! d) solve the velocity field
         CALL csrvec(-1D0,ixBGS%GuvM1T,xzp,yuv)
         ! yuv = buv - Guv*xtilp - GuvM1T*xzp
         CALL inisol(cgtype,mgmres,maxnits,loctolabs,loctolred)
         Prc => ixBGS%PAuv
         xuv = 0
         IF (locgmres) THEN ! solve Auv iteratively with preconditioned gmres
            CALL prcgmres(mgmres,nuv,ixBGS%Auv,Prc,Prc%mlp%perm,xuv,yuv)
         ELSE               ! apply preconditioner of Auv once
            CALL applprc(Prc,xuv,yuv)
         END IF
      ELSEIF (type.eq.2) THEN
         ! c) and d) solve big Saddle Point Problem of
         !    depth averaged pressure and velocity field using a
         !    Simpler preconditioner
         yuvzp = 0
         buvzp(1:nuv) = yuv
         buvzp(nuv+1:nuv+nzp) = bzp
         IF (locgmres) THEN
            CALL inisol(cgtype,mgmres,maxnits,spploctolabs,spploctolred)
            ! CALL inisol(cgtype,mgmres,maxnits,1D-6,1D-6)
            CALL bigsppfgmres(mgmres,nuv+nzp,ixBGS,yuvzp,buvzp)
         ELSE
            yuvzp = buvzp
            CALL apply_bigsppprec (ixbgs, yuvzp, 'SR')
         END IF
         xuv = yuvzp(1:nuv)
         xzp = yuvzp(nuv+1:nuv+nzp)
         ! correct pressure for constant field
!         xzp = xzp - dot_product(xzp,ixBGS%svp3)*ixBGS%svp3
!         xzp = xzp - dot_product(xzp,ixBGS%svp4)*ixBGS%svp4         
      END IF
      xp(nw+1:np) = xzp
      
      ! d) solve the vertical velocity-field
      yw = bp(1:nw)
      CALL csrvec(-1D0,ixBGS%Duv1,xuv,yw)
      CALL lsolve(ixBGS%Aw,yw,xw)
      
      ! e and f) solve the temperature and salinity equation
      yTS = bTS
      CALL csrvec(-1D0,ixBGS%BTSuv,xuv,yTS)
      CALL csrvec(-1D0,ixBGS%BTSw ,xw ,yTS)
!      CALL applprc(ixBGS%PATS,xTS,yTS)
      IF (locgmres) THEN ! solve ATS iteratively with preconditioned gmres
         Prc => ixBGS%PATS
         CALL inisol(cgtype,mgmres,mgmres,loctolabs,loctolred)
         CALL prcfgmres(mgmres,nTS,ixBGS%ATS,Prc,Prc%mlp%perm,xTS,yTS)
      ELSE               ! apply preconditioner of ATS once
         CALL applprc(ixBGS%PATS,xTS,yTS)
      END IF

      ! g) solve the transformed pressure xp
      yw = bw
      xp(1:nw) = 0
      CALL csrvec(-1D0,ixBGS%BwTS,xTS,yw)
      IF (ilu) THEN
         CALL SCfgmres(mgmres,nw,ixBGS,xp(1:nw),yw,3,type)
      ELSE
         CALL usolve(ixBGS%Ap,yw,xp(1:nw))
      END IF

      IF (ilu) THEN
         ! compute correction (in fact application of inverse of U)
         yw = 0; yzp = 0; yTS = 0; yuv = 0; yuv2 = 0;
         ! h) compute correction for xuv and xzp
         CALL csrvec(1D0,ixBGS%Guv1,xp(1:nw),yuv2)  ! yuv2 = Guv1*xp
         IF (type.eq.1) THEN
            CALL applprc(ixBGS%PAuv,yuv,yuv2)       ! yuv = PAuv\yuv2
            yzp = 0
         ELSEIF (type.eq.2) THEN
            buvzp = 0; buvzp(1:nuv) = yuv2;
            !         CALL inisol(cgtype,mgmres,maxnits,spploctolabs,spploctolred)
            !         CALL inisol(cgtype,mgmres,maxnits,1D-6,1D-6)
            !         CALL bigsppfgmres(mgmres,nuv+nzp,ixBGS,yuvzp,buvzp) ! solve bigSPP
            yuvzp = buvzp
            CALL apply_bigsppprec (ixbgs, yuvzp, 'SR')
            yuv = yuvzp(1:nuv)                      ! yuv = SPP\[yuv2;0]
            yzp = yuvzp(nuv+1:nuv+nzp)
            ! correct pressure for constant field
!            yzp = yzp - dot_product(yzp,ixBGS%svp3)*ixBGS%svp3
!            yzp = yzp - dot_product(yzp,ixBGS%svp4)*ixBGS%svp4
         END IF
         ! i) compute correction for xw
         yw = 0; yw2 = 0;
         CALL csrvec(-1D0,ixBGS%Duv1,yuv,yw2)
         CALL lsolve(ixBGS%Aw,yw2,yw);              ! yw = -Aw\(Duv1*yuv)
         ! j) compute correction for xTS
         yTS2 = 0
         CALL csrvec(-1D0,ixBGS%BTSuv,yuv,yTS2)
         CALL csrvec(-1D0,ixBGS%BTSw, yw, yTS2)
         CALL applprc(ixBGS%PATS,yTS,yTS2)          ! yTS = ATS\(-BTSw*yw-BTSuv*yuv)
         ! k) add corrections to x
         xuv = xuv - yuv;
         xp(nw+1:np) = xp(nw+1:np) - yzp;
         xw  = xw  - yw;
         xTS = xTS - yTS;
      END IF
      ! l) compute transformation of the pressure
      xzp = xp(nw+1:np); xp(nw+1:np) = 0
      CALL csrvec(1D0,ixBGS%M1T,xzp,xp);
      xp = xp - dot_product(xp,ixBGS%svp1)*ixBGS%svp1
      xp = xp - dot_product(xp,ixBGS%svp2)*ixBGS%svp2         

      ! finally reorder solution and put it in b
      b(ord) = x

      DEALLOCATE(x)
      DEALLOCATE(yp,yw,yw2,yzp,bzp,bzuv)
      DEALLOCATE(yzuv,yuv,yuv2,yTS,yTS2,yz)
      DEALLOCATE(xtilp,xzp,xzuv,xz)
      DEALLOCATE(yuvzp,buvzp)
      CALL inisol(cgtype,mgmres,maxnits,tolred,tolabs)
 !     CALL iniglb(2)
      END SUBROUTINE apply_bilu3

!**************************************************************
      SUBROUTINE delete_BGSprec(ixBGS,type)
!     frees the memory claimed to store the BGS preconditioner
      USE m_bgsprec
      USE m_wfree
      USE m_sparsekit
      USE m_bgspars
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (bgsprec), POINTER :: ixBGS
      INTEGER :: type ! 1 = depth averaged spp, 2 = modified simpler
!     LOCAL
      INTEGER :: ier

      !  subblocks:
      CALL csrfree(ixBGS%Auv)
      CALL csrfree(ixBGS%Guv)
      CALL csrfree(ixBGS%BuvTS)
      CALL csrfree(ixBGS%Gw )
      CALL csrfree(ixBGS%BwTS)
      CALL csrfree(ixBGS%Duv)
      CALL csrfree(ixBGS%Dw )
      CALL csrfree(ixBGS%BTSuv)
      CALL csrfree(ixBGS%BTSw)
      CALL csrfree(ixBGS%ATS)
      !  dummy matrices
      CALL csrfree(ixBGS%Adumw)
      CALL csrfree(ixBGS%Adump)
      CALL csrfree(ixBGS%Bpdumw)
      CALL csrfree(ixBGS%Bdumwp)
      CALL csrfree(ixBGS%Buvdumw)
      CALL csrfree(ixBGS%Bdumwuv)
      !  depth-averaging
      CALL csrfree(ixBGS%M1)
      CALL csrfree(ixBGS%M2)
      CALL csrfree(ixBGS%M1T)
      CALL csrfree(ixBGS%GuvM1T)
      CALL csrfree(ixBGS%M2Duv)
      !  depth-averaged spp approach
      IF (type.eq.1) THEN
         CALL csrfree(ixBGS%Mzuv)
         CALL csrfree(ixBGS%MzuvT)
         CALL csrfree(ixBGS%MAuv)
         CALL csrfree(ixBGS%MGuv)
         CALL csrfree(ixBGS%MDuv)
         CALL csrfree(ixBGS%MAuvGradDiv)
         CALL csrfree(ixBGS%MCp)
      ELSEIF (type.eq.2) THEN
      !  modified simpler approach
         CALL csrfree(ixBGS%Cp)
      END IF
      !  blocks of Ap
      CALL csrfree(ixBGS%Ap)
      CALL csrfree(ixBGS%InvAp)
      CALL csrfree(ixBGS%Guv1)
      !  block of Aw
      CALL csrfree(ixBGS%Aw)
      CALL csrfree(ixBGS%InvAw)
      CALL csrfree(ixBGS%Duv1)
      !  diagonal blocks of Auv
      IF ((xctsc).or.(type.eq.2)) CALL csrfree(ixBGS%DAuv)
!
!      CALL csrfree(ixBGS%DTS)
!      CALL csrfree(ixBGS%DTSinv)
!      CALL csrfree(ixBGS%DTST)
!      CALL csrfree(ixBGS%DTSinvT)
!
      !  preconditioners
      CALL prcfree(ixBGS%PAuv)
      IF ((xctsc).or.(type.eq.2)) CALL csrfree(ixBGS%PDAuv)
      IF (type.eq.2) CALL prcfree(ixBGS%PCp)
      CALL prcfree(ixBGS%PATS)
      !  preconditioners (2)
      IF (type.eq.1) THEN
         CALL csrfree(ixBGS%PDMAuv)
         CALL prcfree(ixBGS%PMAuv)
         CALL prcfree(ixBGS%PMCp)
      END IF
      !  singular vectors of the pressure
      DEALLOCATE(ixBGS%svp1)
      DEALLOCATE(ixBGS%svp2)
      DEALLOCATE(ixBGS%svp3)
      DEALLOCATE(ixBGS%svp4)
!
      DEALLOCATE(ixBGS)
      NULLIFY(ixBGS)
!
      END SUBROUTINE delete_BGSprec

!**************************************************************
      SUBROUTINE set_prec_parameters(prec,spptype,em,euv,euvGD,ep,eTS,maxiter,rest)
      USE m_bgspars
      USE m_mrilupars
      IMPLICIT none
!     IMPORT/EXPORT
      INTEGER :: prec ! 1 = mrilu, 2 = bilu
      INTEGER :: spptype ! 1 = depth averaged spp, 2 = modified simpler
      REAL    :: em,euv,euvGD,ep,eTS
      INTEGER :: maxiter,rest
      !
      IF (prec.eq.1) THEN
         epsw = em
         mgmres1  = maxiter
         maxnits1 = rest*maxiter
      ELSEIF ((prec.eq.2).AND.(spptype.eq.1)) THEN
         MAuv_epsw = euvGD
         Auv_epsw  = euv
         ATS_epsw  = eTS
         mgmres  = maxiter
         maxnits = rest*maxiter
      ELSEIF ((prec.eq.2).AND.(spptype.eq.2)) THEN
         MAuv_epsw = ep
         Auv_epsw  = euv
         ATS_epsw  = eTS
         mgmres  = maxiter
         maxnits = rest*maxiter
      ELSE
         STOP "unknown prec type in set_prec_parameters"
      END IF
      
      END SUBROUTINE set_prec_parameters      

!**************************************************************

      END MODULE m_bgskit
