!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define diamatrix  anymatrix

#endif

MODULE m_wrtbmtde

CONTAINS

SUBROUTINE wrtbmtde (basenm, a, ad, ier)

USE m_fstrlen
USE m_build
USE m_wrtbbd
USE m_wrtmt

CHARACTER (LEN=*)			, INTENT(IN)		:: basenm
TYPE (csrmatrix)                        , POINTER		:: a
TYPE (diamatrix)                        , POINTER		:: ad
INTEGER					, INTENT(OUT)           :: ier

!     Writes the sparse matrix  A  into the 2 binary files
!        '<basenm>.diag'  and  '<basenm>.offd'
!     The matrix is stored in a Block Diagonal extracted CSR format.
!     The inverse of the block diagonal of A is stored in
!        [ ad%co ]
!     and is written into  '<basenm>.diag'.
!     The off-diagonal elements of A are stored in
!        [A%N, a%beg, a%jco, a%co]
!     and are written into  '<basenm>.offd'.

!     File structure
!     --------------
!     The order of the numbers in the file are
!        1   * integer:  order of matrix [ = A%N ]
!        A%N+1 * integer:  indices indicating beginning of a row
!                        of  A + inv(D)
!        NNZ * integer:  column numbers of nonzero of  A + inv(D)
!                        N.B. Not necessarily in increasing order!!!
!                        {NNZ = a%beg(A%N+1) - 1 + A%N}
!        NNZ * double precision: nonzero values of matrix  A + inv(D)

!     Arguments:
!     ==========
!     basenm   i   Base of the filenames
!     A%N        i   The order of the matrix A
!     Ad%blksiz   i   Number of rows/columns in a block of the block-
!                  diagonal matrix.
!     a%beg     i   a%beg(i) = index in array a%jco or a%co for the first
!                  nonzero element A(i,j) in row i  (1 <= i <= A%N+1)
!                  a%beg(A%N+1) = index of last nonzero element + 1
!     a%jco     i   a%jco(v) = column number of v-th nonzero element
!     a%co      i   a%co(v)  = value of v-th nonzero element
!                  (1 <= v <= a%beg(A%N+1) - 1
!     ad%co     i   ad%co(nz)  = value of i-th diagonal element  inv(D)
!     ier      o   The error code,
!                  =  0   No error
!                  = -3   Error during closing, opening or rewinding
!                         file <basenm>.diag  or  <basenm>.offd.
!                  = -5   Error during writing <basenm>.diag  or
!                         <basenm>.offd.
!                  =-12   No logical unit available to write on
!                         <basenm>.diag  or  <basenm>.offd.

!#enddoc

!     1996-07-17  ,  Doeke de Vries
!     1996-09-24  ,  Mark A. Aves
!     1997-02-27  ,  Doeke D. de Vries  (No Automatic arrays!!!)
!     2003-03-05  ,  Last update (Doeke de Vries)

!-----------------------------------------------------------------------

!     Local parameter:
!     ================

INTEGER, PARAMETER :: maxfnmlen = 80

!     Local variables:
!     ================
INTEGER :: baslen

CHARACTER (LEN=maxfnmlen) :: filenm

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'wrtbmtde'

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

baslen =  MIN (fstrlen(basenm), maxfnmlen-4)

filenm(1:baslen+5) = basenm(1:baslen) // '.diag'

#ifdef DEBUG
PRINT '( A, A )' , 'Block diagonal with:  vsm ', filenm(1:baslen+5)
#endif
CALL wrtbbd (filenm(1:baslen+5), ad, ier)
IF ( ier /= 0 )  GO TO 1000

filenm(1:baslen+5) = basenm(1:baslen) // '.diag'

#ifdef DEBUG
PRINT '( A, A )' , 'Off-diagonals  with:  vsm ', filenm(1:baslen+5)
#endif
CALL wrtmt (filenm(1:baslen+5), .true., a, 2, ier)
IF ( ier /= 0 )  GO TO 1000

!     Normal return:
RETURN

!     Error return:
1000 CONTINUE

!     End of  wrtbmtde
END SUBROUTINE wrtbmtde

END MODULE