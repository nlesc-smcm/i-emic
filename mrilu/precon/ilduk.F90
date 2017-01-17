!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define partmatrix  anymatrix

#endif

#ifdef MKL

#define norm(x) nrm2(x)

#else

#define norm(x) SQRT(SUM(x**2))

#endif

MODULE m_ilduk

CONTAINS
 
SUBROUTINE ilduk (maxLev, compfctr, SingLU, LUTol, a, Part)

USE m_dump
USE m_glbpars
USE m_build
USE m_csrresize

#ifdef MKL
USE mkl95_precision, ONLY: wp => dp
USE mkl95_blas, ONLY: nrm2
#endif

INTEGER					, INTENT(IN)		:: maxLev
DOUBLE PRECISION			, INTENT(IN)            :: compfctr
LOGICAL					, INTENT(IN)            :: SingLU
DOUBLE PRECISION			, INTENT(IN)            :: LUTol
TYPE (csrmatrix)			, POINTER		:: a
TYPE (partmatrix)			, POINTER		:: Part

!     Computes an Incomplete LDU factorization of the Part%N x Part%N matrix  A
!     with maximum Level of fill in, K = MaxLev. (ILU(K))
!     A zero row in the matrix A is not allowed!

!     After the subroutine call the matrices  A, L, D, U, R  satisfy:
!        A = L D U + R
!     Where:
!     L  is a lower triangular matrix, with diagonal elements equal 1.
!     D  is a non-singular diagonal matrix.
!     U  is an upper triangular matrix, with diagonal elements equal 1.
!     R  is the Rest matrix.

!     Arguments:
!     ==========
!     MaxLev   i   The maximum Level of fill in, MaxLev >= 0.
!     Part%N        i   Number of rows/columns in the matrices A, L, D and U.
!     CompFctr i   0<=CompFctr<=1.  Compensation Factor for the diagonal.
!                  The sum of the discarded (dropped) elements in a row
!                  of LD and DU, multiplied by 'CompFctr', are added to
!                  the diagonal element of D.
!     SingLU   i   Singular DU factor allowed in the complete or
!                  incomplete LDU factorization of last block.
!                  Only the last diagonal element of  DU  of the complete
!                  LDU factorization may be "zero"!
!     LUTol    i   Tolerance to determine the singularity of factor U,
!                  in LU-factorisation of last block.
!                  "U is singular" <==> MIN(i::ABS(U(i,i))) <= LUTol * MAX(i::ABS(U(i,i)))
!     A%beg     i   A%beg(i): index in 'A%jco' and 'A%co' of the first
!                  off-diagonal element in row i of matrix A.
!     A%jco     i   A%jco(nz): column number of off-diagonal element A%co(nz).
!     A%co      i   A%co(nz): value of off-diagonal element.
!     Part%offd%beg    o   Part%offd%beg(i): index in 'Part%offd%jco' and 'Part%offd%co' of the first
!                  off-diagonal element in row i of matrix L + (U P).
!     Part%offd%jco    o   Part%offd%jco(nz): column number of off-diagonal element
!                  in row i of matrix L + (U P).
!     Part%offd%co     o   Part%offd%co(nz): value of off-diagonal element of L + (U P).
!     Part%dia     o   The diagonal of  inverse(D).
!     Part%piv    o   The Part%piv vector which characterizes the
!                  permutation matrix  P:  (P x)(k) = x(Part%piv(k)). Here the identity matrix  I.
!     Part%Lnzl     o   Part%Lnzl(i): index in 'Part%offd%jco' and 'Part%offd%co' of the last
!                  non-zero off-diagonal element of matrix L in row 'i'.

!#enddoc

!     Array arguments used as local variables:
!     ========================================
!     ColNr    ColNr(nz), 1 <= nz <= LnzL < i, column number of nz-th
!              non-zero element that should be stored in  L.
!              ColNr(nz), i <= nz <= LnzU <=Part%N, column number of nz-th
!              non-zero element that should be stored in  DU.
!     StorCol  StorCol(j), 1 <= j <= Part%N,
!              = .TRUE.  Element in column 'j' of actual row should be
!                     stored in new L+U.
!              = .FALSE.   Element in column 'j' of actual row should not be
!                     stored in new L+U.
!     Levi     Levi(j), 1 <= j <= Part%N, Level of fill in of (i,j)-element
!              with actual row 'i' and column 'j' of new L+U.
!     Level    Level(nz), 1 <= nz <= Part%offd%nnz, Level of fill in of the
!              non-zero element in L+U.
!     LUi      LUi(j), 1 <= j <= Part%N, value to be stored in column 'j'
!              of actual row of new L+U.

!     Local Parameters:
!     =================

CHARACTER (LEN=*), PARAMETER :: rounam = 'ilduk'

!     Local Variables:
!     ================
!     actrow           Row number of actual row;
!                      in the row-loop:  actrow = row
!     invDi            Value of 1 / D(row,row).
!     LnzL             Last non-zero in L-part of LUi (1:LnzL)
!     LnzU             Last non-zero in U-part of LUi (row:LnzU)
!     MaxDia           MAX(k: 1 <= k <= row: ABS(D(k,k)))
!     MinDia           MIN(k: 1 <= k <= row: ABS(D(k,k)))
!                      = ABS(D(MinRow,MinRow))
!     MinRow           Column number with minimal diagonal element.
!     col              Column number of actual column.
!     Thennz              Number of non-zeros to be stored in actual row
!                      of L or U.
!     Part%offd%nnz            Number of non-zeros in off-diagonal part of L+U.
!     normtmp          norm of A(row,:):
!                         SQRT(SUM(j: 1<=j<=Part%N: A(row,j)**2))
!     nzRow            Index in 'ColNr' of non-zero in actual row 'LUi'.

INTEGER 					:: actrow, i, j, N, minrow, row, col, last, Lnzl, Lnzu, thennz, nzj, nzrow, LUnnz, maxnnz, ier
DOUBLE PRECISION 				:: dropsum, invdi, LUij, normtmp
DOUBLE PRECISION 				:: MaxDia, maxval, MinDia
LOGICAL, ALLOCATABLE, DIMENSION(:)		:: storcol
INTEGER, ALLOCATABLE, DIMENSION(:)		:: colnr, Levi
INTEGER, ALLOCATABLE, DIMENSION(:)		:: Level
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: LUi

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

N = Part%N

ALLOCATE( storcol(1:N), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( colnr(1:N), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( LUi(1:N), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
IF ( maxlev > 0 ) THEN
  ALLOCATE( Levi(1:N), STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
  ALLOCATE( Level(1:5*maxlev*A%nnz), STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
END IF 

!     Clear column search pointer:

storcol = .false.

!     Initialise the permutation vector  Part%piv:
FORALL (i=1:N) Part%piv(i)  = i

LUnnz            = 0
Part%offd%beg(1) = 1

DO  row = 1, N
!        Rows  1,2,...,row-1  have been stored in  Part%OFFD.
  
!        Save row number of actual row (now available outside loop):
  actrow = row
  
!        Copy non-zeros of A(row,1:N) into vector 'LUi':
  
  storcol(row) = .true.
  colnr(row)   = row
  LUi(row)     = 0.0D0
  IF ( maxlev > 0 )  Levi(row)    = 0
  
  Lnzl = 0
  Lnzu = row
  
  normtmp = norm( A%co(A%beg(row):A%beg(row+1)-1) )
  
  DO i = A%beg(row), A%beg(row+1)-1
    col          = A%jco(i)
    storcol(col) = .true.
    LUi(col)     = A%co(i)
    IF ( maxlev > 0 ) Levi(col)    = 0
    IF (col < row) THEN
      Lnzl        = Lnzl + 1
      colnr(Lnzl) = col
    ELSE IF (col > row) THEN
      Lnzu        = Lnzu + 1
      colnr(Lnzu) = col
    END IF
  END DO
  
  IF (normtmp == 0) THEN
    PRINT '(/, A, X, A, /, A, I10, A)' , 'Fatal error in', rounam, '!', 'Row', row, 'of input matrix is zero.'
    CALL dump(__FILE__,__LINE__,'Fatal error in, Row of input matrix is zero.')
  END IF
  
  
!        Compute the vector LUi and then  L(1:row-1) and U(row+1:N):
  
  dropsum = 0
  
  nzrow = 0
  thennz   = 0
  DO WHILE (nzrow < Lnzl)
    nzrow = nzrow + 1
    
!           Look for the minimum column number 'j' in the remaining
!           non-zero entries in 'LUi(1:row-1)':

    nzj=MINLOC(colnr(nzrow:Lnzl),1)+nzrow-1

!   {  j = ColNr(nzj) = MIN(ColNr(nzRow:LnzL))  }
    
!   Swap  ColNr(nzj) <-->  ColNr(nzRow):

    j=colnr(nzj)
    colnr(nzj)   = colnr(nzrow)
    colnr(nzrow) = j

!   {  j = ColNr(nzRow)  }
    
    storcol(j) = .false.
    LUij       = LUi(j)
    
  IF ( maxlev > 0 ) THEN
    IF (Levi(j) <= maxLev) THEN

!     The Level of fill in of L(row,j) is sufficient small.
      
!     Update vector LUi(j+1:N):

      DO i = Part%Lnzl(j)+1, Part%offd%beg(j+1)-1
        
!       Column number of element of U:

        col = Part%offd%jco(i)
        
        IF (storcol(col)) THEN

!         Existing element, add contribution:

          LUi(col)  = LUi(col) - LUij * Part%offd%co(i)
          Levi(col) = MIN(Levi(col), Levi(j)+Level(i)+1)
        ELSE

!         New non-zero, insert into row 'LUi':

          storcol(col) = .true.
          LUi(col)     = - (LUij * Part%offd%co(i))
          Levi(col)    = Levi(j)+Level(i)+1
          IF (col < row) THEN
            Lnzl        = Lnzl + 1
            colnr(Lnzl) = col
          ELSE IF (col > row) THEN
            Lnzu        = Lnzu + 1
            colnr(Lnzu) = col
          END IF
        END IF
      END DO
      thennz        = thennz + 1
      colnr(thennz) = j
    ELSE
      dropsum = dropsum + LUij
    END IF

  ELSE

!   Update vector LUi(j+1:N):

    DO i = Part%lnzl(j)+1, Part%offd%beg(j+1)-1
      
!     Column number of element of U:
  
      col = Part%offd%jco(i)
      
      IF (storcol(col)) THEN

!       Existing element, add contribution:
   
        LUi(col) = LUi(col) - LUij * Part%offd%co(i)

      ELSE

        dropsum = dropsum - LUij * Part%offd%co(i)

      END IF

    END DO

    thennz        = thennz + 1
    colnr(thennz) = j

  END IF

END DO

  
  
!        {  j >= row  }
  
!        Store the non-zero off-diagonal elements of L(row,1:row-1) in
!        matrix LU:
  
  Part%offd%beg(row+1) = LUnnz+1
  IF (LUnnz + thennz > UBOUND(Part%offd%jco,1)) THEN
    Part%offd%nnz=LUnnz
    maxnnz=(LUnnz + thennz)*2+1000
!   100% more space plus an additional offset for small matrices    
!    PRINT '(A, X, A)' , 'Entry 1:', rounam
    CALL csrresize(Part%offd,maxnnz)
  END IF

  Part%offd%jco(LUnnz+1:LUnnz+thennz) = colnr(1:thennz)
  Part%offd%co (LUnnz+1:LUnnz+thennz) = LUi (colnr(1:thennz)) * Part%diA%com(1,colnr(1:thennz))
  IF ( maxlev > 0 )  Level        (LUnnz+1:LUnnz+thennz) = Levi(colnr(1:thennz))
  LUnnz = LUnnz + thennz
  
!        Save index in 'Part%offd%jco' and 'Part%offd%co' of last non_zero off-diagonal
!        in L(row,1:row-1):
  Part%Lnzl(row) = LUnnz
  
  
  if ( maxlev > 0 ) THEN

!   Determine the non-zero off-diagonal elements of  DU(row,row+1:N):
  
    thennz = 0
    DO i = row+1, Lnzu
      col          = colnr(i)
      storcol(col) = .false.
    
      IF (Levi(col) <= maxLev) THEN
!     The Level of fill in of U(row,col) is sufficient small.
        thennz          = thennz + 1
        colnr(row+thennz) = col
      ELSE
        dropsum = dropsum + LUi(col)
      END IF
    END DO
  END IF
  
!        Update new diagonal element and store value into Diag:
  
  storcol(row) = .false.
  LUi(row)     = LUi(row) + dropsum * compfctr


  maxval     = ABS(LUi(row))
  
  IF (row == 1) THEN
    MaxDia = maxval
    MinDia = maxval
    minrow = row
  ELSE IF (maxval > MaxDia) THEN
    MaxDia = maxval
  ELSE IF (maxval < MinDia) THEN
    MinDia = maxval
    minrow = row
  END IF
  
  IF (maxval <= neglgbl) THEN
    IF (SingLU  .ANd. row == N) THEN
!              Issue warning, change last diagonal element and return:
!     (Nearly) Singular matrix DU is allowed, actrow = N,
!     Issue warning, change last diagonal element and return:

      IF (N == 1) MaxDia = 1.0D0
        IF (OutLev >= 1) THEN
          PRINT '(/, A, X, A, A, 3(/, 3X, A), 1P, E12.5, /)', &
               'Warning from', rounam, '!', & 
                'After LDU-factorisation of last submatrix:',  &
               'Last diagonal element is (nearly) zero!', & 
               'Value of this element is changed into', MaxDia
        END IF
      Part%diA%com(1,N) = 1.0D0 / MaxDia

    ELSE
!              .NOT. SingLU .OR. row < N;  issue error message and return
!     (Nearly) Singular matrix is not allowed and nearly zero diagonal
!     element found:
      PRINT '(/, A, X, A, A, /, 3X, A)', 'Error detected in', rounam, '!', 'During (I)LDU-factorisation of last submatrix:'
      IF (actrow /= N) THEN
        PRINT '(3X, I11, A, /)', actrow, '-th diagonal element is (nearly) zero!'
      ELSE
        PRINT '(3X, A, /)', 'Last diagonal element is (nearly) zero!'
      END IF
      CALL dump(__FILE__,__LINE__, '(Nearly) Singular matrix is not allowed and nearly zero diagonal element found')
    END IF
  END IF
  
  invdi   = 1.0D0 / LUi(row)
  Part%diA%com(1,row) = invdi
  
! If necessary, enlarge the matrix Part%OFFD:

  IF ( maxlev > 0 ) THEN
    last = thennz
  ELSE
    last = lnzu-row 
  END IF

  Part%offd%beg(row+1) = LUnnz+1
  IF (LUnnz + last> UBOUND(Part%offd%jco,1)) THEN
    Part%offd%nnz=LUnnz
!   100% more space plus an additional offset for small matrices    
    maxnnz=(LUnnz + last)*2+1000
!    PRINT '(A, X, A)' , 'Entry 2:', rounam
    CALL csrresize(Part%offd,maxnnz)
  END IF

!  Store the non-zero elements of U(row+1:N) into matrix Part%OFFD:

  storcol(colnr(row+1:row+last))    = .false.
  Part%offd%jco(LUnnz+1:LUnnz+last) =      colnr(row+1:row+last)
  Part%offd%co (LUnnz+1:LUnnz+last) = LUi (colnr(row+1:row+last)) * invdi
  IF ( maxlev > 0 ) Level        (LUnnz+1:LUnnz+last) = Levi(colnr(row+1:row+last))

  LUnnz        = LUnnz + last

  Part%offd%beg(row+1) = LUnnz + 1
  
END DO

DEALLOCATE( storcol, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( colnr, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
DEALLOCATE( LUi, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
IF ( maxlev > 0 ) THEN
  DEALLOCATE( Levi, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
  DEALLOCATE( Level, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
END IF 

IF ( MinDia <= LUTol * MaxDia) THEN
!        Matrix  A  is (nearly) singular:
  actrow = minrow
  IF (SingLU  .ANd. actrow == N) THEN
!     Issue warning, change last diagonal element and return:
!     (Nearly) Singular matrix DU is allowed, actrow = N,
    IF (N == 1) MaxDia = 1.0D0
      IF (OutLev >= 1) THEN
        PRINT '(/, A, X, A, A, 3(/, 3X, A), 1P, E12.5, /)', &
        'Warning from', rounam, '!', &
        'After LDU-factorisation of last submatrix:',  &
        'Last diagonal element is (nearly) zero!', &
        'Value of this element is changed into', MaxDia
      END IF
    Part%diA%com(1,N) = 1.0D0 / MaxDia
  ELSE
!              .NOT. SingLU .OR. actrow < N;  issue error message and return
!     (Nearly) Singular matrix is not allowed and nearly zero diagonal
!     element found:
      PRINT '(/, A, X, A, A, /, 3X, A)', 'Error detected in', rounam, '!', 'During (I)LDU-factorisation of last submatrix:'
      IF (actrow /= N) THEN
        PRINT '(3X, I11, A, /)', actrow, '-th diagonal element is (nearly) zero!'
      ELSE
        PRINT '(3X, A, /)', 'Last diagonal element is (nearly) zero!'
      END IF
      CALL dump(__FILE__,__LINE__, '(Nearly) Singular matrix is not allowed and nearly zero diagonal element found')
  END IF
END IF

Part%offd%nnz      = LUnnz
Part%offd%beg(N+1) = LUnnz + 1

END SUBROUTINE ilduk

END MODULE
