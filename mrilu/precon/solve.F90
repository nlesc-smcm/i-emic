!#begindoc
 
#ifndef WITH_UNION

#define prcmatrix  anymatrix
#define partmatrix  anymatrix

#endif

#ifndef MKL

#define dot(x,y) DOT_PRODUCT(x,y)

#endif

MODULE m_solve

CONTAINS

SUBROUTINE solve (P, y)

USE m_dump
USE m_build
USE m_csrvec
USE m_mterrmsg
USE m_diavec
USE m_cscvec

#ifdef WITH_ATLAS 

EXTERNAL :: dgetrs

#else

USE m_dgetrs

#endif

TYPE (prcmatrix)			, POINTER		:: P
DOUBLE PRECISION, DIMENSION(1:P%n)	, INTENT(IN OUT)        :: y

!     This routine solves the multilevel partitioned system
!        P x = y,
!     of dimension P%N*P%N, and returns the solution in  y.

!     Arguments:
!     ==========
!     P%N       i   Number of equations
!     P      	i   the descriptor for the multilevel
!                   partitioned matrix  P.
!     y        	io  In:  The right hand side of the linear system.
!                   Out: The solution of the linear system.

!#enddoc

!     The structure of an MLP matrix is described in the include file
!     build.F90

!     This routine performs first a forward elimination step for a block
!     lower triangular matrix system where the block factors are stored
!     in CSC format.

!     Described below is the algorithm applied to a block lower
!     triangular matrix when the equations are permuted into place.
!     The routine can similarly be applied to a block partitioned
!     upper triangular matrix.

!     Consider the block partitioned (P%N*P%N)-matrix system  L x = b

!        (  I  | 0  | 0   )  ( x1 )     ( b1 )
!        ( ----+----+---- )  ( -- )     ( -- )
!        (  A2 | I  | 0   )  ( x2 )  =  ( b2 )
!        ( ----+----+---- )  ( -- )     ( -- )
!        (  A3 | B3 | C3  )  ( x3 )     ( b3 )

!     Forward eliminating the first block partition gives

!        ( x1 ) = ( b1 )

!        ( I  | 0  ) ( x2 )   ( b2 )   ( A2 )
!        ( ---+--- ) ( -- ) = ( -- ) - ( -- ) ( x1 )
!        ( B3 | C3 ) ( x3 )   ( b3 )   ( A3 )

!     whence updating the right-hand side vector for the next block
!     to solve becomes

!                           ( b2 )   ( A2 ) ( x1 )        ( b2 )
!                         = ( -- ) - ( -- )          ->   ( -- )
!                           ( b3 )   ( A3 )               ( b3 )

!     The system to be solved is now:

!        ( x2 ) = ( b2 )

!        ( C3 ) ( x3 ) = ( b3 ) - ( B3 ) ( x2)   ->   ( b3 )


!     The linear system corresponding with the last partition, here
!     with sub-matrix C3, is solved directly using the L and U factors
!     stored in  C3.


!     Then a backward elimination step for a block upper triangular
!     matrix system is performed, where the block factors are stored
!     in CSR format.

!     Consider the block partitioned (N*N)-matrix system  L x = b

!        (  I | F2 | F3  )  ( x1 )     ( b1 )
!        ( ---+----+---- )  ( -- )     ( -- )
!        (  0 |  I | E3  )  ( x2 )  =  ( b2 )
!        ( ---+----+---- )  ( -- )     ( -- )
!        (  0 |  0 | I   )  ( x3 )     ( b3 )

!     The following backwards elimination of the block partition
!     gives the solution  x:

!        ( x3 )  =  ( b3 )

!        ( x2 )  =  ( b2 ) - ( E3 ) ( x3 )

!        ( x1 )  =  (b1 ) - ( F2 | F3 ) ( x2 )
!                                       ( -- )
!                                       ( x3 )

!     Local Variables:
!     ================
INTEGER 					:: poffs, psz, ier
TYPE (partmatrix), POINTER			:: Par
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: x

CHARACTER (LEN=*), PARAMETER :: rounam = 'solve'

DOUBLE PRECISION, PARAMETER :: one = 1.0D0
DOUBLE PRECISION, PARAMETER :: zero = 0.0D0

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

ALLOCATE( x(1:P%n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

!     Solves the partitioned lower triangular system:
!       (L*(1/D)) z = y
!     where   P =: L D U
!     The solution is returned in  y.

!     Actual partition starts with the first one:
Par => P%mlp%first

!     Loop over linked-list to consider each partition, bar the last:
DO WHILE (ASSOCIATED( Par ))
  
!   Partition offset:
    poffs = Par%off

    SELECT CASE (Par%typ)

    CASE (pldutp)

!     Partition is of type PLDU


!     Extract  Block Size, Partition Size  and the value part of the
!     (block-)diagonal matrix and check if storage type is DIAtp:

!     x(POffs+1:POffs+Par%dia%n) = inv(D)(1:Par%dia%n,1:Par%dia%n) y(POffs+1:POffs+Par%dia%n)

      CALL diavec(.false., one, Par%dia, y(poffs+1:POffs+Par%dia%n), x(poffs+1:POffs+Par%dia%n))

!     Lower triangular part of partition:

      IF (Par%ltr%typ == csctp) THEN

!       Lower triangular block stored in CSC format
  
!       Get the matrix pointers, then solve.
!       The solution is returned in  y, the partition size in  Par%dia%n.

!       y(POffs+Par%ltr%n+1:N) = y(POffs+Par%ltr%n+1:N)
!                              - PL(POffs+Par%ltr%n+1:N,1:Par%ltr%n) x(POffs+1:POffs+Par%ltr%n)

        CALL cscvec (-1.0D0, Par%ltr, x(poffs+1:poffs+Par%ltr%n), y)
      ELSE
!       Symmetric, not implemented, or an illegal matrix/partition type
        CALL dump(__FILE__,__LINE__,'Symmetric, not implemented, or an illegal matrix/partition type')
      END IF

    CASE (pffptp)

#     ifdef DEBUG
        PRINT '(A)' , 'Found Partition with Full Factors!'
#     endif

!     Partition Size:

      psz = Par%fm%n

!     Solves the system  L U x = Y  and returns the solution in Y:

#     ifdef WITH_ATLAS
        CALL dgetrs ('Notranspose', psz, 1, Par%fm%com, psz, Par%piv, y(poffs+1:poffs+psz), psz, ier)
#     else
        CALL dgetrs (psz, Par%fm%com, Par%piv, y(poffs+1:poffs+psz))
#     endif

    CASE (psfptp)

#     ifdef DEBUG
        PRINT '(A)' , 'Found Partition with Sparse Factors!'   
#     endif

!     Solve last diagonal block:

!     The solution is returned in Y.
    
      CALL solldu (Par, y(poffs+1:Par%dia%n))

    CASE DEFAULT
      CALL dump(__FILE__,__LINE__,'Illegal matrix type')
    END SELECT



   Par => Par%next
END DO

! ### SOLVE THE BLOCK UPPER TRIANGULAR SYSTEM  (D U) Z = Y
! ### THE SOLUTION IS RETURNED IN Y.

! Last partition:

Par => P%mlp%last%prev


! *** Loop backwards through linked-list to consider each partition,
!     bar the last
DO WHILE ( ASSOCIATED( Par) )

#ifdef DEBUG
!        Ensure the partition is of type PLDU:
  IF (Par%typ /= pldutp) CALL dump(__FILE__,__LINE__,'Wrong matrix type, not type PLDU')
#endif

!        Partition offset
poffs = Par%off

!        Upper Triangular sub-matrix:

IF (Par%utr%typ == csrtp) THEN
!           Upper Triangular sub-matrix stored in CSR format:
  
!           Get the matrix pointers, then solve
!           The solution is returned in Y, the partition size in  Psz.
! y(POffs+1:POffs+Par%utr%n) = y(POffs+1:POffs+Par%utr%n)
!                         - PU(1:Par%utr%n,POffs+Par%utr%n+1:N) y(POffs+Par%utr%n+1:N)
  CALL csrvec (-1.0D0, Par%utr, y, y(poffs+1:poffs+Par%utr%n))
ELSE
! Symmetric, not implemented, or an illegal matrix/partition type!
  CALL dump(__FILE__,__LINE__,'Symmetric, not implemented, or an illegal matrix/partition type')
END IF

! Extract Block Size, Partition Size and the value part of the
! (block-)diagonal matrix and check if storage type is DIAtp:

! x(POffs+1:POffs+Par%dia%n) = inv(D)(1:Par%dia%n,1:Par%dia%n) y(POffs+1:POffs+Par%dia%n)

  CALL diavec (.false., one, Par%dia, y(poffs+1:POffs+Par%dia%n), x(poffs+1:POffs+Par%dia%n))

! Copy solution X(POffs+1:POffs+Par%dia%n) into Y(POffs+1:POffs+Par%dia%n):

  y(poffs+1:poffs+Par%dia%n) = x(poffs+1:poffs+Par%dia%n)

  Par => Par%prev
END DO

DEALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

END SUBROUTINE solve

!=======================================================================

SUBROUTINE solldu (lu, b)

USE m_dump
USE m_build

#ifdef MKL
USE mkl95_precision, ONLY: wp => dp
USE mkl95_blas, ONLY: dot
#endif

TYPE (partmatrix)			, POINTER		:: LU
DOUBLE PRECISION, DIMENSION(LU%dia%n)	, INTENT(IN OUT)	:: b


!     This routine solves the linear system

!         L 1/D (U P) x = b

!     in which
!     L  is a sparse lower triangular matrix with ones on the diagonal,
!     U  is a sparse upper triangular matrix with ones on the diagonal,
!     D  is a nonsingular diagonal matrix, and
!     P  is a permutation matrix.

!     ??-??-?? , Auke van der Ploeg    Original version
!     96-04-03 , Doeke de Vries        Cleanup

!     Arguments:
!     ==========
!     N        		i  Order of the matrices  D,  L  and  U
!     LU%dia%com    	i  LU%dia%com(i) = value of  1/D(i,i)
!     LU%offd%beg    	i  LU%offd%beg(i) = index in array LU%offd%jco or LU%offd%co for the first
!                 	   nonzero element L(i,j) or U(i,j) in row i  (1<=i<=N)
!                 	   LU%offd%beg(N+1) = index of last nonzero element + 1
!     LU%offd%beg    	o  LU%offd%beg(i): index in array 'LU%offd%jco' and 'LU%offd%co' of the
!                 	   first non-zero off-diagonal element in row i of matrix  L + (U P),  (1<=i<=N).
!                 	   LU%offd%beg(N+1) = index of last nonzero element + 1
!     LU%offd%jco    	i  LU%offd%jco(Nz): LU%offd%comn number of Nz-th non-zero off-diagonal element in matrix  L + (U P).
!                 	   LU%offd%beg(i) <= Nz <= LU%lnzl(i):
!                    	   LU%offd%jco(Nz)  is new LU%offd%comn number
!                 	   LU%lnzl(i)+1 <= Nz <= LU%offd%beg(i+1)-1:
!                    	   LU%offd%jco(Nz)  is LU%offd%comn number before LU%piving.
!     LU%lnzl   	i  LU%lnzl(i) = index in 'LU%offd%jco' and 'LU%offd%co' of the last
!                 	   non-zero off-diagonal element in i-th row of  L.
!     LU%piv    	i  The LU%piv vector which characterizes the permutation matrix  P:  (P x)(k) = x(LU%piv(k)).
!     LU%offd%co     	i  LU%offd%co(Nz)  = value of Nz-th non-zero element (1 <= Nz <= LU%offd%beg(N+1) - 1
!     b        		i  Right hand side of the linear system.
!              		o  Solution vector

!-----------------------------------------------------------------------

!     Local variables:
!     ================
!     i           Row number

INTEGER 					:: i, n, p, q, ier
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)	:: x

#ifdef DEBUG
CHARACTER (LEN=*), PARAMETER :: rounam = 'solLDU'

!     TRACE INFORMATION
PRINT '(A)' , rounam
#endif

n = LU%dia%n

ALLOCATE( x(1:n), STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')

! Solve  x  from  L x = b, store the solution in 'x':

DO i = 1, n
  p = LU%offd%beg(i)
  q = LU%lnzl(i)
  x(i) = b(i) - dot( LU%offd%co(p:q), x(LU%offd%jco(p:q)))
END DO

! Solve  y  from  D (U P) y = x,  or from  (U P) y = (1/D) x, and store the solution  y  in  'b':

DO i = n, 1, -1
! Permute the solution into place
  p = LU%lnzl(i)+1
  q = LU%offd%beg(i+1)-1
  b(LU%piv(i)) = LU%dia%com(1,i) * x(i) - dot( LU%offd%co(p:q), b(LU%offd%jco(p:q)) )
  
END DO

DEALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

END SUBROUTINE solldu

END MODULE