!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define cscmatrix  anymatrix
#define csrdmatrix anymatrix
#define scdematrix anymatrix
#define scbmmatrix anymatrix
#define partmatrix anymatrix
#define mlpmatrix  anymatrix

#endif

MODULE m_wcompr

CONTAINS

RECURSIVE SUBROUTINE wcompr (x)

USE m_dump
USE m_build
USE m_csrresize
USE m_cscresize

TYPE (anymatrix)	, POINTER                     :: x

!     Compress the storage occupied by the descriptor, indicated by 'x', with
!     all the referenced segments from this descriptor.

!     Only implemented for the matrix types:
!     CSC, CSR, DIAtp, FMtp, CSRD, SCBM, PLDU, PFFP and PSFP.

!     Arguments:
!     ==========
!     x      io  In:  Location of matrix descriptor.
!                  Out: New location of the matrix descriptor.

!#enddoc

!     Local Parameters:                                 :
!     =================

CHARACTER (LEN=*), PARAMETER :: rounam = 'wcompr'

!     Local Variables:
!     ================

INTEGER				:: ier
INTEGER 			:: ityp, n
INTEGER, DIMENSION(:), POINTER  :: iarray
TYPE (cscmatrix), POINTER       :: xcsc
TYPE (csrmatrix), POINTER       :: xcsr
TYPE (csrdmatrix), POINTER      :: xcsrd
TYPE (mlpmatrix), POINTER      	:: xmlp
TYPE (partmatrix), POINTER      :: xpart
TYPE (scbmmatrix), POINTER      :: xscbm
TYPE (scdematrix), POINTER      :: xscde

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)', 'Entry:', rounam
#endif

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

ityp = x%typ
  
IF (ityp == csctp) THEN
  xcsc => anytocsc(x)
  n = xcsc%n
  IF (n < 0) CALL dump(__FILE__,__LINE__,'Internal error')
    
  IF (UBOUND(xcsc%beg,1)/=n+1) THEN
    ALLOCATE( iarray(1:n+1), STAT=ier)
    IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
    iarray=xcsc%beg(1:n+1)
    DEALLOCATE( xcsc%beg, STAT=ier)
    IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
    xcsc%beg=>iarray
  END IF
    
  CALL cscresize( xcsc, csrnnz(xcsr) )

ELSE IF (ityp == csrtp) THEN
  xcsr => anytocsr(x)
  n = xcsr%n
  IF (n < 0) CALL dump(__FILE__,__LINE__,'Internal error')
      
  IF (UBOUND(xcsr%beg,1)/=n+1) THEN
    ALLOCATE( iarray(1:n+1), STAT=ier)
    IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
    iarray=xcsr%beg(1:n+1)
    DEALLOCATE( xcsr%beg, STAT=ier)
    IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
    xcsr%beg=>iarray
  END IF
    
  CALL csrresize( xcsr, xcsr%nnz )

ELSE IF (ityp == csrdtp) THEN
  xcsrd => anytocsrd(x)

  n = xcsrd%n
  IF (n < 0) CALL dump(__FILE__,__LINE__,'Internal error')
        
  IF (ASSOCIATED(xcsrd%lotr)) THEN
    IF (UBOUND(xcsrd%lotr,1)/=n) THEN
      ALLOCATE( iarray(1:n), STAT=ier)
      IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
      iarray=xcsrd%lotr(1:n)
      DEALLOCATE( xcsrd%lotr, STAT=ier)
      IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
      xcsrd%lotr=>iarray
    END IF
  END IF
    
  CALL wcompr( csrtoany(xcsrd%offd) ) 
          
ELSE IF (ityp == scdetp) THEN
  xscde => anytoscde(x)

  n = xscde%n
  IF (n < 0) CALL dump(__FILE__,__LINE__,'Internal error')
        
  CALL wcompr( csrtoany(xscde%offd) ) 
          
ELSE IF (ityp == diatp) THEN

ELSE IF (ityp == fmtp) THEN
          
ELSE IF (ityp == mlptp) THEN
  
  xmlp => anytomlp(x)
  xpart => xmlp%first
  DO WHILE (ASSOCIATED( xpart ) )
    CALL wcompr(parttoany(xpart))
    xpart => xpart%next
  END DO
          
ELSE IF (ityp == pffptp) THEN
  xpart => anytopart(x)          

  CALL wcompr(fmtoany(xpart%fm) )

ELSE IF (ityp == pldutp) THEN
  xpart => anytopart(x)          
          
  CALL wcompr(csctoany(xpart%ltr) )
  CALL wcompr(csrtoany(xpart%utr) )
  CALL wcompr(diatoany(xpart%dia) )

ELSE IF (ityp == psfptp) THEN
  xpart => anytopart(x)          
          
  CALL wcompr(csrtoany(xpart%offd) )
  CALL wcompr(diatoany(xpart%dia) )

ELSE IF (ityp == scbmtp) THEN
  xscbm => anytoscbm(x)          

  CALL wcompr(csrtoany(xscbm%a22) )
  CALL wcompr(csrtoany(xscbm%a21) )
  CALL wcompr(csrtoany(xscbm%a12) )
  CALL wcompr(diatoany(xscbm%a11d) )

ELSE

  CALL dump(__FILE__,__LINE__,'Unknown type tag')

END IF

!     Normal Return:
RETURN
      
END SUBROUTINE wcompr

END MODULE