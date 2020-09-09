      MODULE m_sparsekit
!
      CONTAINS

!**************************************************************
      SUBROUTINE invert_ordering(ord,iord)
!     inverts the ordering 'ord' in 'iord'
      IMPLICIT none
!     IMPORT/EXPORT
      INTEGER, DIMENSION(:),INTENT(IN)  :: ord
      INTEGER, DIMENSION(:),INTENT(OUT) :: iord
!     LOCAL
      INTEGER i         ! dummy variable
!     
      DO i = 1,size(ord)
         iord(ord(i))=i
      ENDDO
      END SUBROUTINE invert_ordering

!**************************************************************
      SUBROUTINE reorder_matrix(ixA,ord,iord)
!     reorder csr matrix A according to 'ord'
      USE m_build
      USE m_wacsr
      USE m_wfree
      IMPLICIT none
      
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER          :: ixA
      INTEGER, DIMENSION(:), INTENT(IN)  :: ord,iord
!     LOCAL
      TYPE (csrmatrix), POINTER  :: ixB
      INTEGER              :: dumint,h,i,j,k,lwb,upb,ier,dim,nnzA
      ! dummy variables, counters, bounds
      REAL                 :: dumreal ! swap variable
      ! extract dimension and csr-structure of ixA
      dim = ixA%n
      nnzA = ixA%beg(dim+1) - 1
      ! check compatibility with ord
      IF (SIZE(ord).NE.dim) THEN
         STOP "in reorder_matrix (m_gsolve) ord has wrong dimension"
      ELSEIF (SIZE(iord).NE.dim) THEN
         STOP "in reorder_matrix (m_gsolve) iord has wrong dimension"
      ENDIF
      ! allocate new sparsematrix
      CALL wacsr(dim,nnzA,ixB)
      ! first reorder rows
      ! construct new beg array
      ixB%beg(1) = 1
      ixB%beg(dim+1) = nnzA+1
      DO i = 1,dim
         j = ord(i)
         ixB%beg(i+1) = ixB%beg(i)+ ixA%beg(j+1)-ixA%beg(j)
      ENDDO
      ! construct new jco and co array
      DO i = 1,dim
         lwb = ixB%beg(i)
         upb = ixB%beg(i+1)-1
         k = ixA%beg(ord(i))
         DO j = lwb,upb
            ixB%jco(j) = ixA%jco(k)
            ixB%co(j)  = ixA%co(k)
            k = k+1
         ENDDO
      ENDDO
      ixA%beg  = ixB%beg
      ixA%jco  = ixB%jco
      ixA%co   = ixB%co
      ! reorder columns
      DO i = 1,nnzA
         ixA%jco(i) = iord(ixB%jco(i))
      ENDDO
      ixA%nnz = nnzA
      CALL sort_csr(ixA)

      CALL csrfree(ixB)
      END SUBROUTINE reorder_matrix

!**************************************************************
      SUBROUTINE copy_csr(ixA,ixB)
      USE m_build
      USE m_wacsr
      ! copy A to B
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER :: ixA,ixB
!     LOCAL
      INTEGER :: n,nnzA,ier
      
      n = ixA%n
      nnzA = ixA%beg(n+1)-1
      CALL wacsr(n,nnzA,ixB)
      ixB%beg = ixA%beg
      ixB%jco = ixA%jco
      ixB%co  = ixA%co
      ixB%nnz = nnzA
      END SUBROUTINE copy_csr

!**************************************************************
      SUBROUTINE merge_csr(ixA,ixB,ixC)
      USE m_build
      USE m_wacsr
      ! build matrix C = [A;B]
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER :: ixA,ixB,ixC
!     LOCAL
      INTEGER :: nA,nB,nnzA,nnzB,ier
      
      nA = ixA%n
      nB = ixB%n
      nnzA = ixA%beg(nA+1)-1
      nnzB = ixB%beg(nB+1)-1
      CALL wacsr(nA+nB,nnzA+nnzB,ixC)
      ixC%beg(1:nA)             = ixA%beg(1:nA)
      ixC%beg(nA+1:nA+nB)       = ixB%beg(1:nB) + nnzA
      ixC%beg(nA+nB+1)          = nnzA + nnzB + 1
      ixC%nnz                   = nnzA + nnzB
      ixC%jco(1:nnzA)           = ixA%jco(1:nnzA)
      ixC%jco(nnzA+1:nnzA+nnzB) = ixB%jco(1:nnzB)
      ixC% co(1:nnzA)           = ixA% co(1:nnzA)
      ixC% co(nnzA+1:nnzA+nnzB) = ixB% co(1:nnzB)
      END SUBROUTINE merge_csr

!**************************************************************
      SUBROUTINE sort_csr(ixA)
!     Reorders the rows of a csr matrix such that the in jco
!     the indices on each row are ordered in ascending order
      USE m_build
      USE m_qsorti
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER :: ixA
!     LOCAL
      INTEGER, DIMENSION(:), POINTER :: jcorow
      REAL,    DIMENSION(:), POINTER :: corow
      INTEGER, DIMENSION(:), ALLOCATABLE :: idx,vec
      INTEGER  :: dim, ier, nrow, begrow, endrow, i

      dim = ixA%n
      DO i = 1,dim
         begrow = ixA%beg(i)
         endrow = ixA%beg(i+1)-1
         nrow = endrow-begrow+1
         ALLOCATE(idx(nrow))
         jcorow => ixA%jco(begrow:endrow)
         corow => ixA%co(begrow:endrow)
!         CALL ISORTQX('a',nrow,jcorow,1,idx)
         CALL QSORTI(idx,nrow,jcorow)
         jcorow = jcorow(idx)
         corow = corow(idx)
         DEALLOCATE(idx)
      ENDDO
      NULLIFY(jcorow)
      NULLIFY(corow)
      END SUBROUTINE sort_csr

!**************************************************************
      LOGICAL FUNCTION is_sorted_csr(ixA)
!     True if the csr matrix ixA has rows that are sorted
      USE m_build
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER :: ixA
!     LOCAL
      INTEGER  :: nA, i, j
      
      is_sorted_csr = .TRUE.
      nA = ixA%n
      loop:DO i = 1,nA
         DO j = ixA%beg(i),ixA%beg(i+1)-2
            IF (ixA%jco(j)>ixA%jco(j+1)) THEN
               is_sorted_csr = .FALSE.
               WRITE(6,"(""A is not rowsorted csr"", I6,I6)") i,j
               EXIT loop
            ENDIF
         ENDDO
      END DO loop
      END FUNCTION is_sorted_csr

!**************************************************************
      LOGICAL FUNCTION is_consistent_csr(ixA)
!     debug function: checks if csr matrix is consistent
      USE m_build
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER :: ixA
!     LOCAL
      INTEGER :: i,j,k, n,nnzA
!
      is_consistent_csr = .TRUE.
      IF (.NOT.associated(ixA)) THEN
         WRITE(6,*) "matrix not associated"
         is_consistent_csr = .FALSE.
      ELSE
         n = ixA%n
         IF (ixA%beg(1).EQ.1) THEN
            checkbeg: DO i = 1,n
               IF (ixA%beg(i+1).LT.ixA%beg(i)) THEN
                  WRITE(6,*) "beg array is corrupt: i, beg(i), beg(i+1)", i, ixA%beg(i), ixA%beg(i+1)
                  is_consistent_csr = .FALSE.
                  EXIT checkbeg
               END IF
            END DO checkbeg
         ELSE
            is_consistent_csr = .FALSE.
         END IF
         nnzA = ixA%beg(n+1)-1
         checkjco: DO i = 1,nnzA
            IF (ixA%jco(i).LT.1) THEN
               WRITE(6,*) "jco array is corrupt: i, jco(i)", i, ixA%jco(i)
               is_consistent_csr = .FALSE.
               EXIT checkjco
            END IF
         END DO checkjco
      END IF
      END FUNCTION is_consistent_csr

!**************************************************************
      LOGICAL FUNCTION orthogonal(ixA,ixB)
!     True if the two matrices A and B are orthogonal (i.e) A*B=0
      USE m_build
      USE m_matvec
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER :: ixA,ixB
!     LOCAL
      INTEGER  :: nA, nB, nC, ier
      REAL, DIMENSION(:),ALLOCATABLE :: x, y, z
      REAL :: nrm
      !      
      nA = ixA%n
      nB = ixB%n
      nC = maxval(ixB%jco)
      ALLOCATE(x(nA),y(nB),z(nC))
      CALL random_number(z)
      z = z/DSQRT(dot_product(z,z))
!      CALL matvec(nB,ixB,z,y,ier)
!      CALL matvec(nA,ixA,y,x,ier)
      CALL matvec(nB,1D0,csrtoany(ixB),z,y)
      CALL matvec(nA,1D0,csrtoany(ixA),y,x)
      nrm = dsqrt(dot_product(x,x))
      IF (nrm.LT.1D-10) THEN
         orthogonal = .TRUE.
      ELSE
         orthogonal = .FALSE.
      END IF
      DEALLOCATE(x,y,z)
      END FUNCTION orthogonal

!**************************************************************
      SUBROUTINE extract_submatrix(ixA,ixSubA,rlwb,rupb,clwb,cupb)
!     Extracts the block between rows rlwb and rupb
!     and columns clwb and cupb
!     form the sparse matrix ixA. The subblock is put in ixSubA
      USE m_build
      USE m_wacsr
      USE m_wfree
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER :: ixA,ixSubA
      INTEGER, INTENT(IN) :: rlwb,rupb,clwb,cupb
!     LOCAL
      TYPE (csrmatrix), POINTER :: ixB
      INTEGER  :: nA, nB, nnzA, nnzB, rowA, rowB, jA, jB, ier,i
      
      nA = ixA%n
      nB = rupb - rlwb + 1
      nnzB = ixA%beg(rupb+1)-ixA%beg(rlwb)
      CALL wacsr(nB,nnzB,ixB)
      ixB%beg(nB+1) = nnzB+1
      jB = ixB%beg(1)
      DO rowA = rlwb,rupb
         rowB = rowA-rlwb+1
         loop: DO jA = ixA%beg(rowA),ixA%beg(rowA+1)-1
            IF (ixA%jco(jA).LT.clwb) THEN
               CONTINUE
            ELSE IF (ixA%jco(jA).GT.cupb) THEN
               !EXIT loop 
               !or in case of not row-ordered csr:
               CONTINUE
            ELSE
               ixB%jco(jB) = ixA%jco(jA)-clwb+1
               ixB%co(jB)  = ixA%co(jA)
               jB = jB+1
            END IF
         END DO loop
         ixB%beg(rowB+1) = jB
      END DO
      nnzB = jB-1
      CALL wacsr(nB,nnzB,ixSubA)
      ixSubA%beg = ixB%beg
      ixSubA%beg(nB+1) = jB
      ixSubA%jco = ixB%jco(1:nnzB)
      ixSubA%co  = ixB%co(1:nnzB)
      ixSubA%nnz = nnzB
      CALL csrfree(ixB)
      END SUBROUTINE extract_submatrix

!**************************************************************
      SUBROUTINE extract_block_diagonal(ixA,ixD,blksz)
        ! extracts the block-diagonal of A
        USE m_build
        USE m_wfree
        USE m_wacsr
        ! IMPORT/EXPORT
        TYPE (csrmatrix), POINTER :: ixA, ixD
        INTEGER                   :: blksz
        ! LOCAL
        INTEGER                   :: n,rowA,jA,jD,shift, colD
        !
        n = ixA%n ! size of the problem
        IF (MOD(n,blksz)/=0) THEN
           STOP 'in extract_block_diagonal (sparsekit) size not a multiple of blocksize'
        END IF
        CALL wacsr(n,n*blksz,ixD)
        DO rowA = 1,n
           jA = ixA%beg(rowA)
           shift = MOD(rowA-1,blksz)
           DO WHILE (ixA%jco(jA)<(rowA-shift))
              jA = jA+1
           END DO
           jD = (rowA-1)*blksz+1
           ixD%beg(rowA) = jD
           DO colD = (rowA-shift),(rowA-shift+blksz-1)
              IF (ixA%jco(jA)==colD) THEN
                 ixD%jco(jD) = ixA%jco(jA)
                 ixD%co (jD) = ixA%co (jA)
                 jA = jA+1
                 jD = jD+1
              ELSE
                 ixD%jco(jD) = colD
                 ixD%co (jD) = 0
                 jD = jD+1
              END IF
           END DO
        END DO
        ixD%nnz = blksz*n
        ixD%beg(n+1) = ixD%nnz + 1
      END SUBROUTINE extract_block_diagonal

!**************************************************************
      SUBROUTINE invert_block_diagonal(ixD,ixDinv,blksz)
        ! extracts the block-diagonal of A
        USE m_build
        USE m_wfree
        USE m_wacsr
        ! IMPORT/EXPORT
        TYPE (csrmatrix), POINTER :: ixD, ixDinv
        INTEGER                   :: blksz
        ! LOCAL
        INTEGER                   :: n,m,jD,ier
        REAL                      :: det
        REAL, POINTER             :: a,b,c,d
        !
        n = ixD%n ! size of the problem
        IF (MOD(n,blksz)/=0) THEN
           STOP 'in invert_block_diagonal (sparsekit) size not a multiple of blocksize'
        ELSEIF (ixD%nnz/=blksz*n) THEN
           STOP 'in invert_block_diagonal (sparsekit) incorrect fill pattern'
        END IF
        CALL wacsr(n,n*blksz,ixDinv)
        ixDinv%beg = ixD%beg
        ixDinv%jco = ixD%jco
        ixDinv%nnz = ixD%nnz
        m = INT(n/blksz)
        IF (blksz==2) THEN
           DO i = 1,m
              jD = blksz*blksz*(i-1)
              a => ixD%co(jD+1)
              b => ixD%co(jD+2)
              c => ixD%co(jD+3)
              d => ixD%co(jD+4)
              det = a*d-b*c
              IF (abs(det)>1e-12) THEN
                 ixDinv%co(jD+1) = d/det
                 ixDinv%co(jD+2) = -b/det
                 ixDinv%co(jD+3) = -c/det
                 ixDinv%co(jD+4) = a/det
              ELSE
                 STOP 'in invert_block_diagonal (sparsekit) singular block'
              END IF
           END DO
        ELSE
           STOP 'in invert_block_diagonal (sparsekit) only blksz=2 implemented'
        END IF
      END SUBROUTINE invert_block_diagonal

!**************************************************************
      SUBROUTINE transpose_csr(ixA,n,m,ixB)
!     transposes an nxm csr matrix A, where B = A'
      USE m_build
      USE m_wacsr
      USE m_qsorti
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER :: ixA, ixB
      INTEGER, INTENT(IN) :: n,m
!     LOCAL
      INTEGER                        :: nA, nnzA, ier, i,j
      INTEGER, DIMENSION(:), POINTER :: ico, idx
      INTEGER, DIMENSION(:), POINTER :: rowidx
      REAL,    DIMENSION(:), POINTER :: rowval
!
      nA = ixA%n
      IF (n.NE.nA) THEN
         STOP "in transpose_csr A has wrong dimension"
      END IF
      nnzA = ixA%beg(nA+1)-1
      CALL wacsr(m,nnzA,ixB)
      ixB%beg(m+1) = nnzA+1
      ixB%nnz = nnzA
      ALLOCATE(ico(nnzA))
      ALLOCATE(idx(nnzA))
      DO i = 1,n
         ico(ixA%beg(i):ixA%beg(i+1)-1) = i
      END DO
!      CALL ISORTQX('a',nnzA,ixA%jco,1,idx)
      CALL QSORTI(idx,nnzA,ixA%jco)
      ixB%jco = ico(idx)
      ixB%co  = ixA%co(idx)
      ! construct new beg-array
      ico = ixA%jco(idx)
      i = 1
      j = 1
      DO WHILE (j.LE.nnzA)
         IF (ico(j).EQ.i) THEN
            ixB%beg(i+1) = ixB%beg(i+1)+1
            j = j+1
         ELSEIF (ico(j).GT.i) THEN
            i = i+1
            ixB%beg(i+1) = ixB%beg(i)
         ELSE
      !      WRITE(6,"("" in transpose_csr wrongly sorted icol array"")")
            STOP "in transpose_csr wrongly sorted icol array"
         END IF
      END DO
      ! possibly the last rows of the matrix are empty
      ixB%beg(i+1:m+1) = j
      DEALLOCATE(idx,ico)
      ! finally make rows sorted
      CALL sort_csr(ixB)
      END SUBROUTINE transpose_csr

!**************************************************************
      SUBROUTINE multadd(ixA,w,ixB,ixC,ixD)
!     computes D = A+w*B*C, where A is a square matrix
      USE m_build
      USE m_wacsr
      USE m_wfree
      USE m_qsorti
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER :: ixA,ixB,ixC,ixD
      REAL                      :: w
!      INTEGER                   :: estnnz
!     LOCAL
      TYPE (csrmatrix), POINTER :: ixE
      LOGICAL, DIMENSION(:), ALLOCATABLE :: storcol
      REAL,    DIMENSION(:), ALLOCATABLE :: valcol
      INTEGER :: nA,nB,nC,ier,mA
      INTEGER :: nnzmax, nnzA, nnzrow, rowC
      INTEGER :: rowA,nzA,nzB,nzC,colA,colB,colC
      REAL    :: facB
      INTEGER, DIMENSION(:), POINTER     :: rowE
      INTEGER, DIMENSION(:), ALLOCATABLE :: idx

      nA = ixA%n
      nB = ixB%n
      nC = ixC%n
!      IF (estnnz=0) THEN
         nnzmax = 10*(ixA%nnz+ixB%nnz+ixC%nnz)
!      END IF
      mA = max(maxval(ixA%jco),maxval(ixC%jco))
      CALL wacsr(nA,nnzmax,ixE)
      ALLOCATE(storcol(mA))
      ALLOCATE(valcol(mA)) 
      storcol = .FALSE.
      nnzA = 0
      DO rowA = 1,nA
         IF (nnzA.GT.nnzmax) THEN
            STOP "in multadd not enough space claimed"
         END IF
         ! put entries of A on the row
         nnzrow = 0
         DO nzA = ixA%beg(rowA),ixA%beg(rowA+1)-1
            nnzrow = nnzrow+1
            colA = ixA%jco(nzA)
            ixE%jco(nnzA+nnzrow) = colA
            storcol(colA) = .TRUE.
            valcol(colA)  = ixA%co(nzA)
         END DO
         ! add entries of w*B*C
         DO nzB = ixB%beg(rowA),ixB%beg(rowA+1)-1     ! nonzeros on row 'rowA' in B
            colB = ixB%jco(nzB)    ! column number of nonzero in B
            facB = ixB%co(nzB)     ! common factor for all nonzeros in C
            DO nzC = ixC%beg(colB), ixC%beg(colB+1)-1 ! nonzeros on row 'colB' in C
               colC = ixC%jco(nzC) ! column number of nonzero in C
               IF (storcol(colC).EQ..FALSE.) THEN  ! a new element, insert into new row
                  nnzrow = nnzrow+1
                  ixE%jco(nnzA+nnzrow) = colC
                  storcol(colC) = .TRUE.
                  valcol(colC) = w*facB*ixC%co(nzC)
               ELSE ! an existing element, add up values
                  valcol(colC) = valcol(colC) + w*facB*ixC%co(nzC)
               END IF
            END DO
         END DO
         ! sort elements in row 'rowA' of E
         rowE => ixE%jco(nnzA+1:nnzA+nnzrow)
         ALLOCATE(idx(nnzrow))
!         CALL ISORTQX('a',nnzrow,rowE,1,idx)
         CALL QSORTI(idx,nnzrow,rowE)
         rowE = rowE(idx)
         ixE%co(nnzA+1:nnzA+nnzrow) = valcol(rowE)
         nnzA = nnzA+nnzrow
         ixE%beg(rowA+1) = nnzA+1
         storcol(rowE) = .FALSE. ! clear storcol for next row
         NULLIFY(rowE)
         DEALLOCATE(idx)
      END DO
      ! finally store D
      CALL wacsr(nA,nnzA,ixD)
      ixD%beg = ixE%beg
      ixD%jco = ixE%jco(1:nnzA)
      ixD%co  = ixE%co(1:nnzA)
      ixD%nnz = nnzA
      CALL csrfree(ixE)
      DEALLOCATE(storcol,valcol)
      IF (nnzA>nnzmax) THEN
         WRITE (6,*) "In multadd nnzmax too small"
      END IF
      END SUBROUTINE multadd

!**************************************************************
      SUBROUTINE mult(ixB,ixC,ixD,estnnz)
!     computes D = B*C, where B and C are csr matrices
      USE m_build
      USE m_wacsr
      USE m_wfree
      USE m_qsorti
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER :: ixB,ixC,ixD
      INTEGER                   :: estnnz
!     LOCAL
      TYPE (csrmatrix), POINTER :: ixE
      LOGICAL, DIMENSION(:), ALLOCATABLE :: storcol
      REAL,    DIMENSION(:), ALLOCATABLE :: valcol
      INTEGER :: nB,nC,mC,ier
      INTEGER :: nnzmax, nnzA, nnzrow, rowC
      INTEGER :: rowB,nzB,nzC,colB,colC
      REAL    :: facB
      INTEGER, DIMENSION(:), POINTER     :: rowE
      INTEGER, DIMENSION(:), ALLOCATABLE :: idx

!      WRITE (6,*) "entering mult"
!      IF (.NOT.is_consistent_csr(ixB)) WRITE(6,*) "matrix B is not csr "
!      IF (.NOT.is_consistent_csr(ixC)) WRITE(6,*) "matrix C is not csr "
      nB = ixB%n
      nC = ixC%n
      mC = maxval(ixC%jco)  ! probably the length of mC
!      nnzmax = 3*(ixB%beg(nB+1)+ixC%beg(nC+1)-2)
      nnzmax = estnnz
      ! NB usefull only in this application because of depth averaging
      CALL wacsr(nB,nnzmax,ixE)
      ALLOCATE(storcol(mC))
      ALLOCATE(valcol(mC)) 
      storcol = .FALSE.
      nnzA = 0
      DO rowB = 1,nB
         IF (nnzA.GT.nnzmax) THEN
            STOP "in mult not enough space claimed"
         END IF
         nnzrow = 0
         ! put entries of B*C on the row
         DO nzB = ixB%beg(rowB),ixB%beg(rowB+1)-1     ! nonzeros on row 'rowB' in B
            colB = ixB%jco(nzB)    ! column number of nonzero in B
            facB = ixB%co(nzB)     ! common factor for all nonzeros in C
            DO nzC = ixC%beg(colB), ixC%beg(colB+1)-1 ! nonzeros on row 'colB' in C
               colC = ixC%jco(nzC) ! column number of nonzero in C
               IF (storcol(colC).EQ..FALSE.) THEN  ! a new element, insert into new row
                  nnzrow = nnzrow+1
                  ixE%jco(nnzA+nnzrow) = colC
                  storcol(colC) = .TRUE.
                  valcol(colC) = facB*ixC%co(nzC)
               ELSE ! an existing element, add up values
                  valcol(colC) = valcol(colC) + facB*ixC%co(nzC)                  
               END IF
            END DO
         END DO
         ! sort elements in row 'rowB' of E
         rowE => ixE%jco(nnzA+1:nnzA+nnzrow)
         ALLOCATE(idx(nnzrow)) 
!         CALL ISORTQX('a',nnzrow,rowE,1,idx)
         CALL QSORTI(idx,nnzrow,rowE)
         rowE = rowE(idx)
         ixE%co(nnzA+1:nnzA+nnzrow) = valcol(rowE)
         nnzA = nnzA+nnzrow
         ixE%beg(rowB+1) = nnzA+1
         storcol(rowE) = .FALSE. ! clear storcol for next row
         NULLIFY(rowE)
         DEALLOCATE(idx)
      END DO
      ! finally store D
      CALL wacsr(nB,nnzA,ixD)
      ixD%beg = ixE%beg
      ixD%jco = ixE%jco(1:nnzA)
      ixD%co  = ixE%co(1:nnzA)
      ixD%nnz = nnzA
      CALL csrfree(ixE)
      DEALLOCATE(storcol,valcol)
      IF (nnzA>nnzmax) THEN
         WRITE (6,*) "In mult nnzmax too small"
      END IF

      END SUBROUTINE mult

!**************************************************************
      SUBROUTINE build_inverse(ixU,ixQ)
!     constructs the exact inverse Q of the upper triangular matrix U
!     such that: Q*U = U*Q = Id
!     this subroutine only works for the very special case that
!     U has at most one off-diagonal on each row
      USE m_build
      USE m_wacsr
      USE m_wfree
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER     :: ixU,ixQ
!     LOCAL
      TYPE (csrmatrix), POINTER     :: ixA
      INTEGER                       :: i,j,k, nnzA
!
      ! estimate of the amount of nonzeros in the inverse
      nnzA = (ixU%n)*((ixU%n)/(ixU%jco(2)-1))
      CALL wacsr(ixU%n,nnzA,ixA)
      ixA%nnz = nnzA
      ixA%beg(ixA%n+1) = ixA%nnz+1
      j = 1   ! count nnz
      DO i = 1, ixU%n
         IF (j>nnzA) THEN
            WRITE(6,*) "error: not enough space in build_inverse (sparsekit)", j, nnzA
            STOP
         END IF
         k = i
         ixA%beg(i) = j
         DO WHILE ((ixU%beg(k+1)-ixU%beg(k)).EQ.2)
            k = ixU%beg(k)
            ixA%jco(j) = ixU%jco(k)
            ixA%co (j) = 1/ixU%co(k)
            j = j+1
            k = ixU%jco(k+1)
         END DO
         k = ixU%beg(k)
         ixA%jco(j) = ixU%jco(k)
         ixA%co (j) = 1/ixU%co(k)
         j = j+1
      END DO
      nnzA = j-1
      ixA%beg(ixU%n+1) = nnzA+1
      CALL wacsr(ixU%n,nnzA,ixQ)
      ixQ%beg = ixA%beg(1:ixU%n+1)
      ixQ%jco = ixA%jco(1:nnzA)
      ixQ% co = ixA% co(1:nnzA)
      ixQ%nnz = nnzA
      CALL csrfree(ixA)
      END SUBROUTINE build_inverse

!**************************************************************
      SUBROUTINE build_singular_matrix(ixU,n,m,ixM)
!     The matrix U of dimension n*m (n<m) is supposed to be
!     an upper two diagonal matrix of full rank.
!     This subroutine constructs the matrix M such that
!     such that U*M' = 0 = M*U'
      USE m_build
      USE m_wacsr
      USE m_wfree
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER     :: ixU, ixM
      INTEGER                       :: n,m
!     LOCAL
      TYPE (csrmatrix), POINTER     :: ixB ! dummy matrix
      INTEGER                       :: ier,i,j,k,svs
      LOGICAL, DIMENSION(:),POINTER :: track
      REAL                          :: nrm,sf
!
      ! constuct the singular vectors of U
      CALL wacsr(n,m,ixB)
      ixB%nnz = n
      ixB%beg(ixB%n+1) = ixB%nnz+1
      ALLOCATE(track(ixU%n))
      track = .FALSE.
      svs = 0 ! count singular vectors
      j = 0   ! count nnz
      ixB%co = 1;
      DO i = 1, ixU%n
         IF (.NOT.track(i)) THEN
            svs = svs+1
            k = i
            sf = 1
            DO WHILE (k.LE.ixU%n)
               j = j+1
               track(k) = .TRUE.
               k = ixU%beg(k)
               ixB%jco(j) = ixU%jco(k)
               ixB%co(j)  = sf
               sf = sf*abs(ixU%co(k))/abs(ixU%co(k+1))
               k = ixU%jco(k+1)
            END DO
            j = j+1
            ixB%jco(j) = k
            ixB%co(j) = sf
            ixB%beg(svs+1) = j+1
         END IF
      END DO
      CALL wacsr(svs,j,ixM)
      ixM%beg = ixB%beg(1:svs+1)
      ixM%jco = ixB%jco
      ixM%co  = ixB%co
      ixM%nnz = j
      DO i = 1,ixM%n
         j = ixM%beg(i)
         k = ixM%beg(i+1)
         nrm = dsqrt(dot_product(ixM%co(j:k-1),ixM%co(j:k-1)))
!	 nrm = 1
!	 nrm = abs(ixM%co(k-1))
         if (nrm>0D0) then
            ixM%co(j:k-1) = ixM%co(j:k-1)/nrm 
         else
            WRITE(*,*) " error in  build_singular_matrix: nrm=0, i,j,k  ", i,j,k
         end if
      END DO
      CALL csrfree(ixB)
      DEALLOCATE(track)

      END SUBROUTINE build_singular_matrix

!**************************************************************
      SUBROUTINE lsolve(ixL,b,x)
!     solves equation Lx = b, where L is lower triangular csr matrix
      USE m_build
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER        :: ixL
      REAL, DIMENSION(1:), INTENT(IN)  :: b
      REAL, DIMENSION(1:), INTENT(OUT) :: x
!     local
      INTEGER :: i,nz,diag,dim,ier
      REAL    :: sum

      dim = ixL%n
      DO i = 1,dim
         diag = ixL%beg(i+1)-1
         IF (ixL%jco(diag).EQ.i) THEN
            sum = b(i)
            DO nz = ixL%beg(i),diag-1
               sum = sum - ixL%co(nz)*x(ixL%jco(nz))
            END DO
            x(i) = sum/ixL%co(diag)
         ELSE
!            CALL print_sparse_matrix(ixL)
            STOP "in lsolve L is not lower triangular or not sorted"
         END IF
      END DO
      
      END SUBROUTINE lsolve

!**************************************************************
      SUBROUTINE dsolve(ixD,b,x)
!     solves equation Dx = b, where D is diagonal csr matrix
      USE m_build
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER              :: ixD
      REAL, DIMENSION(1:), INTENT(IN)  :: b
      REAL, DIMENSION(1:), INTENT(OUT) :: x
!     local
      INTEGER :: i,nz,diag,dim,ier
      REAL    :: sum

      dim = ixD%n
      DO i = 1,dim
         IF ((ixD%beg(i).EQ.i).AND.(ixD%jco(i).EQ.i)) THEN
            x(i) = b(i)/ixD%co(i)
         ELSE
!            CALL print_sparse_matrix(ixD)
            STOP "in dsolve D is not diagonal"
         END IF
      END DO
      END SUBROUTINE dsolve

!**************************************************************
      SUBROUTINE usolve(ixU,b,x)
!     solves equation Ux = b, where U is upper triangular csr matrix
      USE m_build
      IMPLICIT none
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER              :: ixU
      REAL, DIMENSION(1:), INTENT(IN)  :: b
      REAL, DIMENSION(1:), INTENT(OUT) :: x
!     local
      INTEGER :: i,nz,diag,dim,ier
      REAL    :: sum

      dim = ixU%n
      DO i = dim,1,-1
         diag = ixU%beg(i)
         IF (ixU%jco(diag).EQ.i) THEN
            sum = b(i)
            DO nz = diag+1,ixU%beg(i+1)-1
               sum = sum - ixU%co(nz)*x(ixU%jco(nz))
            END DO
            x(i) = sum/ixU%co(diag)
         ELSE
            STOP "in usolve U is not upper triangular or not sorted"
         END IF
      END DO
      END SUBROUTINE usolve

!**************************************************************
      SUBROUTINE print_sparse_matrix(ixA)
      USE m_build
!     IMPORT/EXPORT
      TYPE (csrmatrix), POINTER          :: ixA
!     LOCAL
      INTEGER                      :: i,j,n,ier
!
      n = ixA%n
      WRITE(6,"("" Writing matrix: "")")
      WRITE(6,"("" nnz, ico, jco, co "")")
      DO i = 1,n
         DO j = ixA%beg(i), ixA%beg(i+1)-1
            WRITE(6,"(I6,I7,I7,"" "",F16.6)") j,i,ixA%jco(j),ixA%co(j)
         END DO
      END DO
      WRITE(6,"("" End matrix "")")
      END SUBROUTINE print_sparse_matrix

!**************************************************************
      SUBROUTINE write_sparse_matrix(ixA,fname)
!     write the matrix A to a file in asci format
      USE m_build
      IMPLICIT none
!     INPUT/OUTPUT
      TYPE (csrmatrix), POINTER :: ixA
      character(len=*)          :: fname     
!     LOCAL
      INTEGER                   :: i,ier, n
!
      n = ixA%n
!      OPEN(10,file = 'A.asc')
!      WRITE(10,*) n
!      DO i = 1, n+1 
!         WRITE(10,*) ixA%beg(i)
!      END DO
!      DO i = 1, ixA%beg(n+1)-1
!         WRITE(10,*) ixA%jco(i)
!      END DO 
!      DO i = 1, ixA%beg(n+1)-1
!         WRITE(10,*) ixA%co(i)
!      END DO
!      CLOSE(10) 
      OPEN(10,file=fname//'.info')
      OPEN(11,file=fname//'.beg')
      OPEN(12,file=fname//'.jco')
      OPEN(13,file=fname//'.co')
      WRITE(10,*) n,ixA%beg(n+1)-1
      DO i = 1,n+1
        write(11,*) ixA%beg(i)
      ENDDO
      DO i = 1,ixA%beg(n+1)-1
        write(12,*) ixA%jco(i)
      ENDDO
      DO i = 1,ixA%beg(n+1)-1
        write(13,*) ixA%co(i)
      ENDDO
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      END SUBROUTINE write_sparse_matrix

!**************************************************************
      SUBROUTINE print_vector(x)
!     IMPORT/EXPORT
      REAL, DIMENSION(:) :: x
!     LOCAL
      INTEGER            :: i,lb,ub

      WRITE(6,"("" Writing vector: "")")
      lb = lbound(x,1)
      ub = ubound(x,1)
      DO i = lb, ub
         WRITE(6,"(I5,"" "",F16.8)") i, x(i)
      END DO
      WRITE(6,"("" End vector "")")
      END SUBROUTINE print_vector

!**************************************************************
      SUBROUTINE mycsrvec(alpha,A,x,y)
        ! computes y => y+alpha*A*x
        USE m_build
        IMPLICIT none
        ! IMPORT/EXPORT
        TYPE (csrmatrix), POINTER :: A
        REAL, DIMENSION(:)        :: x
        REAL, DIMENSION(1:A%n)    :: y
        REAL                      :: alpha
        ! LOCAL
        INTEGER :: r,nz
        REAL    :: sum
        !
        write(6,*) "entering mycsrvec"
        write(6,*) "A%n = ", A%n, ", A%nnz = ", A%nnz
        DO r = 1,A%n
           IF ((A%n==34560).AND.(r>32507).AND.(r<32515)) THEN
              WRITE(6,*) " in mycsrvec: r = ", r, ", beg(r) = ", A%beg(r), ", beg(r+1) =", A%beg(r+1)
           END IF
           sum = 0D0
           DO nz = A%beg(r),A%beg(r+1)-1
              IF ((A%n==34560).AND.(r>32507).AND.(r<32515)) THEN
                 WRITE(6,*) "nz =", nz, ", co(nz) =",A%co(nz),", jco(nz) =",A%jco(nz)
                 WRITE(6,*) "x  =", x(A%jco(nz))
              END IF
              sum = sum + A%co(nz)*x(A%jco(nz))
           END DO
           y(r) = y(r) + alpha*sum
        END DO
        write(6,*) "mycsrvec finished"
      END SUBROUTINE mycsrvec

!**************************************************************
      END MODULE m_sparsekit
