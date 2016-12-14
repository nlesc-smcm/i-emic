!#begindoc

#ifndef WITH_UNION

#define csrmatrix  anymatrix
#define partmatrix  anymatrix
#define diamatrix  anymatrix
#define scdematrix  anymatrix

#endif

MODULE m_csslpar

CONTAINS
 
SUBROUTINE csslpar (NEqDon, NEqNotDon, BlkSiz, S, Part)

USE m_dump
USE m_build
USE m_wacsr
USE m_wapsfp
USE m_wcompr
USE m_wfree
USE m_glbpars
USE m_prcpars
USE m_vispars
USE m_wrtbldu
USE m_ioerrmsg
USE m_ilduk
USE m_incldup
USE m_scd2csr

INTEGER			, INTENT(IN)		:: NEqDon
INTEGER			, INTENT(IN)            :: NEqNotDon
INTEGER			, INTENT(IN)            :: BlkSiz
TYPE (scdematrix)	, POINTER		:: S
TYPE (partmatrix)	, POINTER		:: Part

!     Computes and Stores Sparse Last Partition into a new area.
!     The storage of the matrix  S  is released.

!     Arguments:
!     ==========
!     NEqDon   i   Number of rows/columns that have been entered in the
!                  Multilevel Partitioned LDU matrix.
!     NEqNotDon    i   Number of rows/columns in the last Schur-complement
!                  matrix  S.
!     BlkSiz   i   Number of rows/columns in a diagonal block.
!     S      io  In:  Location of the Schur-complement matrix  S.
!                  Out: --- undefined.
!     Part   o   Location of descriptor of the newly
!                  allocated last partition with the L, D and U factors.
!#enddoc

!     Local Parameters:
!     =================

CHARACTER (LEN=*), PARAMETER :: rounam = 'csslpar'

!     Local variables:
!     ================

INTEGER 					:: ier
DOUBLE PRECISION 				:: begtim, cumtim, endtim
TYPE (csrmatrix), POINTER           		:: B

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Get the locations of the constituent parts of the
!     Schur Complement, Block Diagonal extracted, matrix:

!     Request workspace for a copy,
!     B, in CSR format of the matrix  S  and return the locations of
!     the constituent parts:
CALL wacsr (NEqNotDon, NEqNotDon*BlkSiz + csrnnz(S%offd), B)

!     Store the the Schur-Complement, Diagonal extracted, matrix S into
!     the matrix B, stored in CSR format:
CALL scd2csr (S, B)

!     Release workspace occupied by matrix S:
CALL scdefree(S)

!     Allocate a PSFP partition [Part]. (Descriptor has been initialised!)

!     This is an estimate of nnz by its upperbound, a better and smaller estimate is possible

    CALL wapsfp (NEqNotDon, INT((NEqNotDon**2)/10), NEqDon, Part)

!     Begin timing:
    CALL CPU_TIME(begtim)

    IF (ilutype == 9) THEN
!     Use ILUT type factorization:
      
!     Make an Incomplete LDU factorization of the last Schur-
!     complement, B, using partial pivoting and store the factors
!     L, D  and  U  into the SCD type matrix of the partition
!     'Part':

      CALL incldup (compfctr, droptol, cpivtol, singlu, lutol, B,  Part)
      
    ELSE
!     Use ILU(K) type factorization, with  K = ILUType = 0,1,2,3,...:
        
!     Make an Incomplete LDU factorization of the last Schur-complement, B, using partial pivoting and store the factors
!     L, D  and  U  into the SCD type matrix of the partition 'Part':

      CALL ilduk (ilutype, compfctr, singlu, lutol, B, Part)
        
    END IF
        
    CALL CPU_TIME(endtim)
    IF (outlev >= 3) THEN
      cumtim = endtim - begtim
      PRINT '(A, X, F7.2, /)', 'Time ILDU factorization of last block:', cumtim
    END IF
        
!   Free workspace for  B:
    CALL csrfree(B)
        
    IF (visildu) THEN
!     Visualisation of Incomplete LDU factorization required.
          
      IF (outlev >= 2) THEN
        PRINT '(A, X, A, 3(/, 3X, A), /)' ,  &
                'Visualize L-, D- and U-factors of last',  &
                'Schur-complement with:', 'vsm  ILDUfactor.D  ,',  &
                'vsm  ILDUfactor.L  and', 'vsm  ILDUfactor.U'
      END IF
          
      CALL wrtbldu ('ILDUfactor', Part, ier)
      IF (ier /= 0) CALL ioerrmsg ('ILDUfactor', 'wrtbldu', rounam, ier, 0)
          
    END IF
        
        
!   Compress the storage by adjusting the required size:
    CALL wcompr (parttoany(Part))
        
!   End of  csslpar
  END SUBROUTINE csslpar

END MODULE