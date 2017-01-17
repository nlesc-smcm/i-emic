!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define cscmatrix  anymatrix
#define csrdmatrix anymatrix
#define scdematrix anymatrix
#define scbmmatrix anymatrix
#define diamatrix  anymatrix

#endif

MODULE m_wennz

CONTAINS

RECURSIVE FUNCTION wennz (x) RESULT(res)

USE m_dump
USE m_build

TYPE (anymatrix), POINTER 	:: x
INTEGER		                :: res

!     Extract the order and the number of nonzeros in the matrix
!     indicated by 'x'.
!     Only implemented for the matrix types:
!     CSC, CSR, CSRD, DIAtp, FMtp  and  SCBM!

!     Arguments:
!     ==========
!     x      i   Location of the Matrix descriptor.
!     Wennz     o   The number of nonzeros in the matrix indicated by
!                  'x'.

!#enddoc

!     Local Parameters:
!     =================

CHARACTER (LEN=*), PARAMETER :: rounam = 'wennz'

!     Local Variables:
!     ================

INTEGER 			:: ityp
TYPE (scbmmatrix), POINTER      :: xscbm
TYPE (scdematrix), POINTER      :: xscde
TYPE (diamatrix), POINTER      	:: xdia
TYPE (cscmatrix), POINTER      	:: xcsc
TYPE (csrmatrix), POINTER      	:: xcsr
TYPE (csrdmatrix), POINTER      :: xcsrd

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)', 'Entry:', rounam
#endif

ityp = x%typ
  
IF (ityp == csctp) THEN
  xcsc => anytocsc(x)
  res   = cscnnz(xcsc)
ELSE IF (ityp == csrtp) THEN
  xcsr => anytocsr(x)
  res   = csrnnz(xcsr)
ELSE IF (ityp == csrdtp) THEN
  xcsrd => anytocsrd(x)
  res = wennz(diatoany(xcsrd%dia)) + wennz(csrtoany(xcsrd%offd))
ELSE IF (ityp == scdetp) THEN
  xscde => anytoscde(x)
  res = wennz(diatoany(xscde%dia)) + wennz(csrtoany(xscde%offd))
ELSE IF (ityp == diatp) THEN
  xdia => anytodia(x)
  res = xdia%blksiz * x%n
ELSE IF (ityp == fmtp) THEN
  res = x%n*x%n
ELSE IF (ityp == scbmtp) THEN
  xscbm => anytoscbm(x)
  res = wennz(diatoany(xscbm%a11d)) + wennz(csrtoany(xscbm%a12)) + wennz(csrtoany(xscbm%a21)) + wennz(csrtoany(xscbm%a22))
ELSE
    
  PRINT '(A, 2X, A, /, 3X, A, I11, 3X, A)' , 'Internal error in', rounam,  &
      'Storage type', ityp, 'not implemented!'
  CALL dump(__FILE__,__LINE__,'Not implementd')
END IF
  
END FUNCTION wennz

END MODULE
