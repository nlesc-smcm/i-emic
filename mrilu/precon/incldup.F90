!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define partmatrix anymatrix

#endif

#ifdef MKL

#define norm(x) nrm2(x)

#else

#define norm(x) SQRT(SUM(x**2))

#endif

MODULE m_incldup

CONTAINS
 
SUBROUTINE incldup (compfctr, DropTol, CPivTol, singLU, LUtol, a, Part)

USE m_dump
USE m_glbpars
USE m_build
USE m_csrresize

#ifdef MKL
USE mkl95_precision, ONLY: wp => dp
USE mkl95_blas, ONLY: nrm2
#endif

DOUBLE PRECISION				, INTENT(IN)    :: compfctr
DOUBLE PRECISION				, INTENT(IN)    :: DropTol
DOUBLE PRECISION				, INTENT(IN)    :: CPivTol
LOGICAL						, INTENT(IN)    :: singLU
DOUBLE PRECISION				, INTENT(IN)    :: LUtol
TYPE (csrmatrix)				, POINTER	:: a
TYPE (partmatrix)				, POINTER	:: Part

!     This subroutine makes a complete or an Incomplete
!     LDU factorization of the A%N x A%N matrix  A  using partial
!     Pivoting  (ILUT).
!     A zero row in the matrix A is not allowed!

!     If DropTol = 0:
!     A complete LDU factorization of A is made.
!     If the last diagonal element, D(A%N,A%N), in the factorization is
!     nearly zero and a singular LU-factorization is allowed, SingLU,
!     then this diagonal element is replaced by
!     MAX(ABS(D(i,i)) : 1 <= i < A%N).

!     If DropTol > 0:
!     An Incomplete LDU factorization of A is made.
!     This factorization is computed in the same (row-oriented) manner
!     as the LDU factorization except after each  row, say 'i', of LD
!     and DU has been calculated, all off-diagonal entries in that row,
!     except the entries at the positions of the nonzero elements of A,
!     which are smaller in  magnitude than the local drop tolerance are
!     "dropped" from the rows  L(i,:) or U(i,:).
!     The local drop tolerance is:
!        DropTol * ||A(i,:)||
!     The diagonal of the upper triangular factor DU, i.e. the pivot
!     element, will never be dropped even if it is too small.  Any zero
!     pivot entry is replaced by the local drop tolerance in an attempt
!     to avoid a singular factor.
!     It is possible to use various diagonal compensation ILDU's such
!     as MILDU with the argument 'CompFctr'.

!     After the subroutine call the matrices  A, L, D, U, R  and the
!     matrix  P  satisfy:
!        A P' = L D U + R
!     <==>
!        A = L D (U P) + R P
!     Where:
!     L  is a lower triangular matrix, with diagonal elements equal 1.
!     D  is a non-singular diagonal matrix.
!     U  is an upper triangular matrix, with diagonal elements equal 1.
!     P  is a permutation matrix, with  inverse(P) = P'.
!     R  is the Rest matrix.

!     The off-diagonal entries in row 'i' of the rest matrix R are less
!     than  DropTol * ||A(i,:)||.
!     Furthermore  R(i,i) = - CompFctr * SUM(j : j~=i : R(i,j))

!     Arguments:
!     ==========
!     A%N        i   Number of rows/columns in the matrices A, L, D and U.
!     CompFctr i   0<=CompFctr<=1.  Compensation Factor for the diagonal.
!                  The sum of the discarded (dropped) elements in a row
!                  of LD and DU, multiplied by 'CompFctr', are added to
!                  the diagonal element of D.
!     DropTol  i   A non-negative scalar used as the drop tolerance for
!                  the incomplete LDU factorization.
!     CPivTol  i   Permutation tolerance. The two columns 'i' and 'j'
!                  are permuted only when:
!                     ABS(DU(i,j)) * CPivTol > ABS(DU(i,i)
!     SingLU   i   Singular DU factor allowed in the complete or
!                  incomplete LDU factorization of last block.
!                  Only the last diagonal element of  DU  of the complete
!                  LDU factorization may be "zero"!
!     LUTol    i   Tolerance to determine the singularity of factor U,
!                  in LU-factorisation of last block.
!                  "U is singular" <==>
!                  MIN(i::ABS(U(i,i))) <= LUTol * MAX(i::ABS(U(i,i)))
!     A%beg     i   A%beg(i): index in 'A%jco' and 'A%co' of the first
!                  off-diagonal element in row i of matrix A.
!     A%jco     i   A%jco(nz): column number of off-diagonal element
!                  A%co(nz).
!     A%co      i   A%co(nz): value of off-diagonal element.
!     Part%offd%beg    o   Part%offd%beg(i): index in 'Part%offd%jco' and 'Part%offd%co' of the first
!                  off-diagonal element in row i of matrix L + (U P).
!     Part%offd%jco    o   Part%offd%jco(nz): column number of off-diagonal element
!                  in row i of matrix L + (U P).
!     Part%offd%co     o   Part%offd%co(nz): value of off-diagonal element of L + (U P).
!     Part%dia    o   The diagonal of  inverse(D).
!     Part%piv    o   The pivot vector which characterizes the
!                  permutation matrix  P:  (P x)(k) = x(Part%piv(k)).
!                  Part%piv(New column number) = Old column number.
!     Part%Lnzl     o   Part%Lnzl(i): index in 'Part%offd%jco' and 'Part%offd%co' of the last
!                  non-zero off-diagonal element of matrix L in row 'i'.

!#enddoc

!     Array arguments used as local variables:
!     ========================================
!     ColNr    ColNr(nz), 1 <= nz <= Lnzl < i, column number of nz-th
!              non-zero element that should be stored in  L.
!              ColNr(nz), i <= nz <= Lnzu <=A%N, column number of nz-th
!              non-zero element that should be stored in  DU.
!     StorCol  StorCol(j), 1 <= j <= A%N,
!              = Drop    Element in column 'j' of actual row should not
!                        be stored in new L+U.
!              = Stornz  Element in column 'j' of actual row should be
!                        stored in new L+U.
!              = Orignz  Element in column 'j' of actual row in the
!                        original matrix A is nonzero and should be
!                        stored in new L+U.
!     LUi      LUi(j), 1 <= j <= A%N, value to be stored in column 'j'
!              of actual row of new L+U.
!     InvPiv   The inverse permutation of the one stored in 'Part%piv'.
!              InvPiv(Old column number) = New column number.


!     Global Parameters:
!     ==================

!     Local Parameters:
!     =================

CHARACTER (LEN=*), PARAMETER :: rounam = 'incldup'


INTEGER, PARAMETER :: drop = 0
INTEGER, PARAMETER :: stornz = 1
INTEGER, PARAMETER :: orignz = 2

!     Local Variables:
!     ================
!     actrow           Row number of actual row;
!                      in the i-loop:  actrow = i
!     invDi            Value of 1 / D(i,i).
!     Lnzl             Last non-zero in L-part of LUi (1:Lnzl)
!     Lnzu             Last non-zero in U-part of LUi (i:Lnzu)
!     MaxDia           MAX(k: 1 <= k <= i: ABS(D(k,k)))
!     MinDia           MIN(k: 1 <= k <= i: ABS(D(k,k)))
!                      = ABS(D(MinRow,MinRow))
!     MinRow           Column number with minimal diagonal element.
!     NewCol           New column number of actual column: InvPiv(A%jco)
!     Thennz              Number of non-zeros to be stored in actual row
!                      of L or U.
!     Lu%nnz            Number of non-zeros in off-diagonal part of L+U.
!     nrm             norm of A(i,:):
!                         SQRT(SUM(j: 1<=j<=A%N: A(i,j)**2))
!     NrzDiag          Number of zeros in diagonal.
!     nzA              Index of non-zero in representation of matrix A.
!     nzLU             Index of non-zero in representation of off-
!                      diagonals of matrix L+U.
!     nzRow            Index in 'ColNr' of non-zero in actual row 'LUi'.

INTEGER 					:: actrow, i, j, jmax, minrow, newcol, Lnzl, lnzu, thennz, nrzdiag, nz, nza, nzj, nzLU, nzmax, nzrow, LUnnz, maxnnz, ier
DOUBLE PRECISION 				:: dropsum, invdi, LUij, nrm
DOUBLE PRECISION 				:: absLUii, absLUij, MaxDia, maxval, MinDia

INTEGER, ALLOCATABLE, DIMENSION(:)		:: storcol, colnr, invpiv
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: LUi

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

ALLOCATE( storcol(1:A%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( colnr(1:A%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( invpiv(1:A%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
ALLOCATE( LUi(1:A%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Clear column search pointer:
storcol = A%n

!     Initialise the permutation vectors  InvPiv  and  Part%piv:
FORALL (i=1:A%n) invpiv(i) = i
Part%piv = invpiv

nrzdiag = 0

LUnnz     = 0
Part%offd%beg(1) = 1

DO  i = 1, A%n
!        Rows  1,2,...,i-1  have been stored in  LU.
  
!        Save row number of actual row (now available outside loop):
  actrow = i
  
!        Copy non-zeros of (A P')(i,1:n) into vector 'LUi':
  
  storcol(i) = orignz
  colnr(i)   = i
  LUi(i)     = 0.0D0
  
  Lnzl = 0
  lnzu = i
  
  DO nza = A%beg(i), A%beg(i+1)-1
    newcol          = invpiv(A%jco(nza))
    storcol(newcol) = orignz
    LUi(newcol)     = A%co(nza)
    IF (newcol < i) THEN
      Lnzl        = Lnzl + 1
      colnr(Lnzl) = newcol
    ELSE IF (newcol > i) THEN
      lnzu        = lnzu + 1
      colnr(lnzu) = newcol
    END IF
  END DO

  nrm = norm( A%co(A%beg(i):A%beg(i+1)-1) )

! The following line was activated for Utrecht computations with extreme
! coefficients in some rows of the matrix in order to avoid that connections
! are dropped whereas corresponding connections in other rows are not
! dropped. In addition to this MRILU should do no row scaling.
!         nrm=1
  
  IF (nrm == 0) THEN
    PRINT '(/, A, 2X, A, /, 3X, A, I8, X, A)' ,  &
        'Fatal error in', rounam, '!', 'Row', i, 'of input matrix is zero.'
    CALL dump(__FILE__,__LINE__,'Fatal error, Row of input matrix is zero.')
  END IF
  
  
!        Compute the vector LUi and then  L(1:i-1) and U(i+1:N):
  
  dropsum = 0
  
  nzrow = 0
  thennz   = 0
  DO WHILE (nzrow < Lnzl)
    nzrow = nzrow + 1
    
!           Look for the minimum column number 'j' in the remaining
!           non-zero entries in 'LUi(1:i-1)':

    nzj=MINLOC(colnr(nzrow:Lnzl),1)+nzrow-1

!   {  j = ColNr(nzj) = MIN(ColNr(nzRow:Lnzl))  }
    
!   Swap  ColNr(nzj) <-->  ColNr(nzRow):

    j            = colnr(nzj)
    colnr(nzj)   = colnr(nzrow)
    colnr(nzrow) = j

!   {  j = ColNr(nzRow)  }
    
    LUij       = LUi(j)
    
    IF ( ABS(LUij) > DropTol * nrm .OR.  storcol(j) == orignz) THEN
!              The element is large enough to store
    
!              Update vector LUi(j+1:N):
      DO nzLU = Part%Lnzl(j)+1, Part%offd%beg(j+1)-1
       
!     Column number of element of U = (U P) P':
      newcol = invpiv(Part%offd%jco(nzLU))
      
      IF (storcol(newcol) /= drop) THEN
!       Existing element, add contribution:
        LUi(newcol) = LUi(newcol) - LUij * Part%offd%co(nzLU)
      ELSE
!       New non-zero, insert into row 'LUi':
        storcol(newcol) = stornz
        LUi(newcol)     = - (LUij * Part%offd%co(nzLU))
        IF (newcol < i) THEN
          Lnzl        = Lnzl + 1
          colnr(Lnzl) = newcol
        ELSE IF (newcol > i) THEN
          lnzu        = lnzu + 1
          colnr(lnzu) = newcol
        END IF
      END IF
    END DO
    thennz        = thennz + 1
    colnr(thennz) = j
  ELSE
    dropsum = dropsum + LUij
  END IF
  storcol(j) = drop
END DO

!        {  j >= i  }

!        Store the non-zero off-diagonal elements of L(i,1:i-1) in
!        matrix Part%offd:

IF (LUnnz + thennz > UBOUND(Part%offd%jco,1)) THEN
  Part%offd%nnz =LUnnz
! 10% more space plus an additional offset for small matrices    
  maxnnz=(LUnnz + thennz)*1.1+1000
  CALL csrresize(Part%offd,maxnnz)
END IF

! Save index in 'Part%offd%jco' and 'Part%offd%co' of last non_zero off-diagonal in L(i,1:i-1):

Part%offd%jco(LUnnz+1:LUnnz+thennz) = colnr(1:thennz)
Part%offd%co (LUnnz+1:LUnnz+thennz) = LUi(colnr(1:thennz)) * Part%dia%com(1,colnr(1:thennz))
LUnnz = LUnnz + thennz
Part%Lnzl(i) = LUnnz


!        Determine the non-zero off-diagonal elements of  DU(i,i+1:N)
!        and the new pivot element in  DU(i,i:N):

absLUii = ABS(LUi(i))
nzmax   = i
maxval  = absLUii

thennz = 0
DO nz = i+1, lnzu 
  j          = colnr(nz)
  absLUij    = ABS(LUi(j))
  
  IF (      absLUij > DropTol * nrm .OR.  storcol(j) == orignz) THEN
!   The element is large enough to store
    thennz          = thennz + 1
    colnr(i+thennz) = j
  
    IF (      (absLUij > maxval) .AND. (absLUij * CPivTol > absLUii) ) THEN
      nzmax  = i + thennz
      maxval = absLUij
    END IF
  ELSE
    dropsum = dropsum + LUi(j)
  END IF
  storcol(j) = drop
END DO


! {  MaxVal = ABS(LUi(jMax)) ,  jMax = ColNr(nzMax)  }
! jMax is the pivot column in row 'i' of  (D U).

! Swap  ColNr(nzMax) <--> ColNr(i):

jmax         = colnr(nzmax)
colnr(nzmax) = colnr(i)
colnr(i)     = jmax


! Update new diagonal element and store value into Diag:

storcol(i) = drop
LUi(jmax)  = LUi(jmax) + dropsum * compfctr
maxval     = ABS(LUi(jmax))

IF (i == 1) THEN
  MaxDia = maxval
  MinDia = maxval
  minrow = i
ELSE IF (maxval > MaxDia) THEN
  MaxDia = maxval
ELSE IF (maxval < MinDia) THEN
  MinDia = maxval
  minrow = i
END IF

IF (maxval <= neglgbl) THEN
    IF ((singLU) .AND. (i == A%n)) THEN
!     Issue warning, change last diagonal element and return:
!     (Nearly) Singular matrix DU is allowed, actrow = N,
      IF (A%n == 1) MaxDia = 1.0D0
        IF (outlev >= 1) THEN
          PRINT '(/, A, X, A, A, 3(/, 3X, A), 1P, E12.5, /)',  &
                'Warning from', rounam, '!',  &
                'After LDU-factorisation of last submatrix:',  &
                'Last diagonal element is (nearly) zero!',  &
                'Value of this element is changed into', MaxDia
      END IF
      Part%dia%com(1,A%n) = 1.0D0 / MaxDia
    ELSE IF ((singLU) .AND. (DropTol > 0.0D0)) THEN

!     DropTol > 0  AND  i < A%N

      LUi(jmax) = (1.0D-4 + DropTol) * nrm
      nrzdiag   = nrzdiag + 1

    ELSE
!     issue error message and return
!     (Nearly) Singular matrix is not allowed and nearly zero diagonal
!     element found:
      PRINT '(/, A, X, A, A, /, 3X, A)', 'Error detected in', rounam, '!',  &
            'During (I)LDU-factorisation of last submatrix:'
      IF (actrow /= A%n) THEN
        PRINT '(3X, I11, A, /)', actrow, '-th diagonal element is (nearly) zero!'
      ELSE
        PRINT '(3X, A, /)', 'Last diagonal element is (nearly) zero!'
      END IF
      CALL dump(__FILE__,__LINE__,'(Nearly) Singular matrix is not allowed and nearly zero diagonal element found')
    END IF
END IF

  invdi   = 1.0D0 / LUi(jmax)
  Part%dia%com(1,i) = invdi


! Store the non-zero elements of (U P)(i+1:N) into matrix Part%offd:

  IF (LUnnz + thennz > UBOUND(Part%offd%jco,1)) THEN
    Part%offd%nnz =LUnnz
!   10% more space plus an additional offset for small matrices    
    maxnnz=(LUnnz + thennz)*1.1+1000
    CALL csrresize(Part%offd,maxnnz)
  END IF

  Part%offd%jco(LUnnz+1:LUnnz+thennz) = Part%piv(colnr(i+1:i+thennz))
  Part%offd%co (LUnnz+1:LUnnz+thennz) = LUi     (colnr(i+1:i+thennz)) * invdi
  LUnnz = LUnnz + thennz

  Part%offd%beg(i+1) = LUnnz + 1

! Update the pivot vector 'Part%piv' and its inverse 'InvPiv' to
! reflect the interchange of column 'i' and the pivot
! column 'jMax':

  IF (jmax /= i) THEN
    j           = Part%piv(jmax)
    Part%piv(jmax) = Part%piv(i)
    Part%piv(i)    = j
  
    invpiv(j)           = i
    invpiv(Part%piv(jmax)) = jmax
  END IF
END DO

IF ( MinDia <= LUtol * MaxDia) THEN
!        Matrix  A  is (nearly) singular:
  actrow = minrow
    IF ((singLU) .AND. (actrow == A%n)) THEN
!     Issue warning, change last diagonal element and return:
!     (Nearly) Singular matrix DU is allowed, actrow = N,
      IF (A%n == 1) MaxDia = 1.0D0

        IF (outlev >= 1) THEN
          PRINT '(/, A, X, A, A, 3(/, 3X, A), 1P, E12.5, /)',  &
                'Warning from', rounam, '!',  &
                'After LDU-factorisation of last submatrix:',  &
                'Last diagonal element is (nearly) zero!',  &
                'Value of this element is changed into', MaxDia
      END IF
      Part%dia%com(1,A%n) = 1.0D0 / MaxDia
    ELSE IF ((singLU) .AND. (DropTol > 0.0D0)) THEN
!           The diagonal element has been changed already.
    ELSE
!     issue error message and return
!     (Nearly) Singular matrix is not allowed and nearly zero diagonal
!     element found:
      PRINT '(/, A, X, A, A, /, 3X, A)', 'Error detected in', rounam, '!',  &
            'During (I)LDU-factorisation of last submatrix:'
      IF (actrow /= A%n) THEN
        PRINT '(3X, I11, A, /)', actrow, '-th diagonal element is (nearly) zero!'
      ELSE
        PRINT '(3X, A, /)', 'Last diagonal element is (nearly) zero!'
      END IF
      CALL dump(__FILE__,__LINE__,'(Nearly) Singular matrix is not allowed and nearly zero diagonal element found')
    END IF
END IF


  IF (nrzdiag > 0) THEN
!        Incomplete LDU factorization:
    IF (outlev >= 1) THEN
      PRINT '(/, A, X, A, A, /, 3X, A)', 'Warning from', rounam, '!',  &
          'During the ILDU-factorisation of last submatrix:'
      IF (nrzdiag == 1) THEN
        PRINT '(3X, A, 1X, A)', '1  diagonal element has been',  &
            'replaced with local drop tolerance!'
      ELSE
        PRINT '(3X, I4, 2X, A, 1X, A)',  &
            nrzdiag, 'diagonal elements have been',  &
            'replaced with local drop tolerance!'
      END IF
      PRINT '(3X, A, /)', 'The preconditioner might not work!'
    END IF
  END IF
  
  
!     Normal return:
  Part%offd%beg(A%n+1) = LUnnz + 1
  Part%offd%nnz =LUnnz

  DEALLOCATE( storcol, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
  DEALLOCATE( colnr, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
  DEALLOCATE( invpiv, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
  DEALLOCATE( LUi, STAT=ier )
  IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')


END SUBROUTINE incldup

END MODULE