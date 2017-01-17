!#begindoc

#ifndef WITH_UNION

#define diamatrix  anymatrix

#endif

MODULE m_wrtbbd

CONTAINS
 
SUBROUTINE wrtbbd (filenm, ad, ier)

USE m_build
USE m_getunit

CHARACTER (LEN=*)		, INTENT(IN)		:: filenm
TYPE (diamatrix)		, POINTER         	:: ad
INTEGER				, INTENT(OUT)           :: ier

!     Writes a Block Diagonal, stored by column, matrix  Ad  in
!     Compressed Sparse Row (CSR) format to the binary, i.e.
!     unformatted, file 'filenm'.

!     The representation of the matrix A is stored, columnwise, in
!     the argument:
!        [ ad%com ]

!     File structure
!     --------------
!     The file 'filenm' will contain 4 records, with
!       record 1   1   * integer:          order of the matrix A (= N).
!       record 2   N+1 * integer:          indices indicating beginning
!                                          of a row of A.
!       record 3   NNZ * integer:          column numbers of nonzeros
!                                          of A.
!                                          N.B. In increasing order.
!                                          {NNZ = begAd(NBlock+1) - 1}
!       record 4   NNZ * double precision: nonzero values of matrix A.

!     Arguments:
!     ==========
!     filenm   i   Filename
!     AD%N        i   The number of rows/colums in the matrix A.
!     Ad%blksiz   i   Number of rows/columns in a diagonal block.
!     ad%com     i   ad%com(v)  = value of v-th element in block Ad(i,i),
!                  ( ( (i-1)*Ad%blksiz**2+1 <= v < i*Ad%blksiz**2 ),
!                    1 <= i <= AD%N / Ad%blksiz ).
!     ier      o   The error code,
!                  =  0   No error
!                  = -3   Error during closing, opening or rewinding
!                         file 'filenm'
!                  = -5   Error during writing 'filenm'
!                  =-12   No logical unit available to write on
!                         'filenm'.

!#enddoc

!     1996-11-18  ,  Mark A. Aves
!     2003-03-05  ,  Last update (Doeke de Vries)


!     Local Variables:
!     ================
INTEGER :: i, j, iblk, nblock, iunit

#ifdef DEBUG

CHARACTER (LEN=*), PARAMETER :: rounam = 'wrtbbd'

!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Get an free logical unit number for writing:

CALL getunit (iunit, ier)
IF ( ier /= 0 )  GO TO 988

OPEN  (iunit, FILE = filenm, STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', ERR = 997)
REWIND (iunit, ERR = 997)

!     Compute number of blocks in the matrix A:
nblock = ad%n / ad%blksiz

!     Write the data to file:
WRITE (iunit, ERR = 995) ad%n
WRITE (iunit, ERR = 995) ( ( (iblk-1)*ad%blksiz**2 + i*ad%blksiz + 1, i = 0, ad%blksiz-1 ), iblk = 1, nblock ),  &
    nblock*ad%blksiz**2 + 1
WRITE (iunit, ERR = 995) ( ( (iblk-1)*ad%blksiz + MOD(i-1,ad%blksiz) + 1, i = 1, ad%blksiz**2 ), iblk = 1, nblock )
WRITE (iunit, ERR = 995) ( ( ( ad%com((iblk-1)*ad%blksiz+i,j), i = 1, ad%blksiz ), j = 1, ad%blksiz ), iblk = 1, nblock )

CLOSE (iunit, ERR = 997)

!     Normal return:
ier = 0
RETURN

!     Error return:
988  ier = ier - 7
995  ier = ier - 2
997  ier = ier - 3

!     End of  wrtbbd
END SUBROUTINE wrtbbd

END MODULE