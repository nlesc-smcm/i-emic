!#begindoc
 
#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define cscmatrix  anymatrix
#define diamatrix  anymatrix
#define fmmatrix   anymatrix
#define csrdmatrix anymatrix
#define scdematrix anymatrix
#define scbmmatrix anymatrix
#define mlpmatrix  anymatrix
#define partmatrix anymatrix
#define pldumatrix anymatrix
#define pffpmatrix anymatrix
#define prcmatrix  anymatrix

#endif

MODULE m_wfree

CONTAINS

RECURSIVE SUBROUTINE anyfree (x)

USE m_dump
USE m_build

TYPE (anymatrix) 	, POINTER                  :: x

!     Releases storage occupied by 'x'

!     Arguments:
!     ==========
!     x   	io  In:  Location of Matrix descriptor.
!                   Out: -- (undefined value)

!#enddoc

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

IF (x%typ == csctp) THEN
 
  CALL cscfree(anytocsc(x))
     
ELSE IF (x%typ == csrtp) THEN

  CALL csrfree(anytocsr(x))
      
ELSE IF (x%typ == diatp) THEN

  CALL diafree(anytodia(x))
          
ELSE IF (x%typ == fmtp) THEN

  CALL fmfree(anytofm(x))
          
ELSE IF (x%typ == csrdtp) THEN

  CALL csrdfree(anytocsrd(x))
     
ELSE IF (x%typ == scdetp) THEN

  CALL scdefree(anytoscde(x))
          
ELSE IF (x%typ == scbmtp) THEN

  CALL scbmfree(anytoscbm(x))
            
ELSE IF (x%typ == mlptp) THEN
  
  CALL mlpfree(anytomlp(x))
            
ELSE IF (x%typ == pldutp) THEN

  CALL pldufree(anytopart(x))

ELSE IF (x%typ == pffptp) THEN

  CALL pffpfree(anytopart(x))
            
ELSE IF (x%typ == psfptp) THEN 

  CALL psfpfree(anytopart(x))

ELSE IF (x%typ == prctp) THEN
  
  CALL prcfree(anytoprc(x))

ELSE

  CALL dump(__FILE__,__LINE__,'Unknown type tag')

END IF
                
END SUBROUTINE anyfree



SUBROUTINE csrfree (x)

USE m_dump
USE m_build

TYPE (csrmatrix) 	, POINTER	:: x

INTEGER					:: ier

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

DEALLOCATE( x%beg, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
  
DEALLOCATE( x%jco, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
     
DEALLOCATE( x%co, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

DEALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

END SUBROUTINE



SUBROUTINE cscfree (x)

USE m_dump
USE m_build

TYPE (cscmatrix) 	, POINTER	:: x

INTEGER					:: ier

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

DEALLOCATE( x%beg, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
  
DEALLOCATE( x%jco, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
     
DEALLOCATE( x%co, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

DEALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
     
END SUBROUTINE


SUBROUTINE diafree (x)

USE m_dump
USE m_build

TYPE (diamatrix) 	, POINTER	:: x

INTEGER					:: ier

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

DEALLOCATE( x%com, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

DEALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
     
END SUBROUTINE



SUBROUTINE fmfree (x)

USE m_dump
USE m_build

TYPE (fmmatrix) 	, POINTER	:: x

INTEGER					:: ier

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

DEALLOCATE( x%com, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

DEALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
     
END SUBROUTINE


SUBROUTINE csrdfree (x)

USE m_dump
USE m_build

TYPE (csrdmatrix) 	, POINTER	:: x

INTEGER					:: ier

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

IF (x%typ /= csrdtp) CALL dump(__FILE__,__LINE__,'Invalid subtype')

CALL csrfree(x%offd)
CALL diafree(x%dia)

DEALLOCATE( x%lotr, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

DEALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
     
END SUBROUTINE



SUBROUTINE scdefree (x)

USE m_dump
USE m_build

TYPE (scdematrix) 	, POINTER	:: x

INTEGER					:: ier

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

IF (x%typ == csrdtp) THEN
! the matrix is actually of the type scde
  CALL csrdfree(scdetocsrd(x))
  RETURN
END IF

IF (x%typ /= scdetp) CALL dump(__FILE__,__LINE__,'Invalid subtype')

CALL csrfree(x%offd)
CALL diafree(x%dia)

DEALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
     
END SUBROUTINE



SUBROUTINE scbmfree (x)

USE m_dump
USE m_build

TYPE (scbmmatrix) 	, POINTER	:: x

INTEGER					:: ier

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

IF (x%typ /= scbmtp) CALL dump(__FILE__,__LINE__,'Invalid subtype')

CALL diafree(x%a11d)
CALL csrfree(x%a12)
CALL csrfree(x%a21)
CALL csrfree(x%a22)
            
DEALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
     
END SUBROUTINE


RECURSIVE SUBROUTINE prcfree (x)

USE m_dump
USE m_build

TYPE (prcmatrix) 	, POINTER	:: x

INTEGER					:: ier

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

IF (x%typ /= prctp) CALL dump(__FILE__,__LINE__,'Invalid subtype')

DEALLOCATE( x%scale, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

DEALLOCATE( x%perrb, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

IF ( ASSOCIATED(x%aro) ) CALL anyfree(x%aro)

IF ( ASSOCIATED(x%mlp) ) CALL mlpfree(x%mlp)

DEALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
     
END SUBROUTINE



RECURSIVE SUBROUTINE mlpfree (x)

USE m_dump
USE m_build

TYPE (mlpmatrix), POINTER	:: x
INTEGER				:: ier
TYPE (partmatrix), POINTER 	:: xdesc, xptr

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

DEALLOCATE( x%perm, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

xptr => x%first
DO WHILE (ASSOCIATED( xptr ) )

  xdesc => xptr
  xptr => xptr%next
  CALL partfree(xdesc)
 
END DO

DEALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
     
END SUBROUTINE


SUBROUTINE partfree (x)

USE m_dump
USE m_build

TYPE (partmatrix), POINTER	:: x

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

IF (x%typ == pldutp) THEN
! the matrix is actually of the type scde
  CALL pldufree(x)
  RETURN
END IF

IF (x%typ == pffptp) THEN
! the matrix is actually of the type scde
  CALL pffpfree(x)
  RETURN
END IF

IF (x%typ == psfptp) THEN
! the matrix is actually of the type scde
  CALL psfpfree(x)
  RETURN
END IF

END SUBROUTINE


SUBROUTINE pldufree (x)

USE m_dump
USE m_build

TYPE (partmatrix) 	, POINTER	:: x

INTEGER					:: ier

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

IF (x%typ /= pldutp) CALL dump(__FILE__,__LINE__,'Invalid subtype')

CALL diafree(x%dia)
CALL cscfree(x%ltr)
CALL csrfree(x%utr)

DEALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
     
END SUBROUTINE


SUBROUTINE pffpfree (x)

USE m_dump
USE m_build

TYPE (partmatrix) 	, POINTER	:: x

INTEGER					:: ier

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

IF (x%typ /= pffptp) CALL dump(__FILE__,__LINE__,'Invalid subtype')

CALL fmfree(x%fm)

DEALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
     
END SUBROUTINE


SUBROUTINE psfpfree (x)

USE m_dump
USE m_build

TYPE (partmatrix) 	, POINTER	:: x

INTEGER					:: ier

IF ( .NOT. ASSOCIATED(x) ) CALL dump(__FILE__,__LINE__,'Unassociated pointer')

IF (x%typ /= psfptp) CALL dump(__FILE__,__LINE__,'Invalid subtype')

CALL diafree(x%dia)
CALL csrfree(x%offd)

DEALLOCATE( x%piv, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
              
DEALLOCATE( x%lnzl, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
                
DEALLOCATE( x, STAT=ier )
IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')

END SUBROUTINE



END MODULE