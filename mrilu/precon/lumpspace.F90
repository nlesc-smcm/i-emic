!#begindoc
 
#ifndef WITH_UNION

#define scdematrix  anymatrix
#define diamatrix  anymatrix

#endif

MODULE m_lumpspace

CONTAINS
 
SUBROUTINE lumpspace (S, CSpace, RSpace)

USE m_dump
USE m_build
USE m_wadia
USE m_wfree
USE m_glbpars
USE m_prcpars
USE m_dgeco
USE m_dgedi

TYPE (scdematrix)				, POINTER		:: S
DOUBLE PRECISION, DIMENSION(1:S%N)		, INTENT(OUT)           :: CSpace
DOUBLE PRECISION, DIMENSION(1:S%N)		, INTENT(OUT)           :: RSpace

!     Calculate the column and row lump space from the Block Diagonal
!     of the Schur Complement, Block Diagonal extracted, matrix  S,
!     referenced via  S.

!     Arguments:
!     ==========
!     S%N    	i   Number of rows/columns in (sub-)matrix  S.
!     S      	i   Location of the descriptor of the Schur-
!              	    Complement, (block-) Diagonal Extracted, matrix  S.
!     CSpace   	o   CSpace(r): "Lump Space tolerance" / MAX(ABS(Sd(r,:))),
!                   with  1<=r<=S%N  and  Sd  is the diagonal part of S.
!     RSpace   	o   RSpace(c): "Lump Space tolerance" / MAX(ABS(Sd(:,c))),
!                   with  1<=c<=S%N  and  Sd  is the diagonal part of S.

!#enddoc

CHARACTER (LEN=*), PARAMETER :: rounam = 'lumpspace'


!     Local Variables:
!     ================
INTEGER 					:: vec, ier
INTEGER, DIMENSION(:), ALLOCATABLE		:: temp
INTEGER 					:: BlkSiz, iBlk, NBlk, nsingBlk
INTEGER 					:: i, BlkSizout
DOUBLE PRECISION 				:: absval, rcond
DOUBLE PRECISION, DIMENSION(:,:), POINTER  	:: Blk
TYPE (diamatrix), POINTER  			:: Bd

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:',  rounam
#endif


!     Initialisation:

CSpace = 0
RSpace = 0
nsingBlk = 0
BlkSiz=S%dia%BlkSiz
  
!     Compute and check the value of  NBlk:
  NBlk   = S%N / BlkSiz
  
  IF (NBlk*BlkSiz /= S%N) THEN
!        Error if 'S%N' and 'NBlk' do not satisfy:
!        "S%N is an integer multiple of BlkSiz"
    PRINT '( /, A, 2X, A, /, 2(3X, A, I10, /) )' ,  &
        'Internal error in', rounam, '!',  &
        'Block size:                           ', BlkSiz,  &
        'inconsistent with number rows/columns:', S%N
    CALL dump(__FILE__,__LINE__,'Internal error')
  END IF
  
  IF (BlkSiz /= 1) THEN
!        Size of diagonal blocks > 1:
    
!        Request Int and DP workspace segments [Pivot, DWork] for the
!        subroutines  dgeco  and  dgedi:
    ALLOCATE( temp(1:BlkSiz), STAT=ier)
    IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
    
!        Request the workspace for a copy of the Block Diagonal matrix:
    CALL wadia (S%N, BlkSiz, Bd)
    
!        Copy the block diagonal matrix, of S, into Bd:
    Bd%com = S%dia%com
    
    vec = 0
    
    DO iBlk = 1, NBlk
!           {  Blk = Bd%com + (iBlk-1)*BlkSiz**2 ,
!              Vec = (iBlk-1)*BlkSiz             }
      
!           Factor the matrix and estimate condition number:
      Blk => Bd%com(1:BlkSiz,(iBlk-1)*BlkSiz+1:iBlk*BlkSiz)

      CALL dgeco (Blk, BlkSiz, temp, rcond )
      
      IF (rcond > 1.0D-8) THEN
        
!              Compute the inverse of a matrix

        CALL dgedi (Blk, BlkSiz, temp)
        
!              Compute Column and Row lump space:
        FORALL (i=1:BlkSiz)
          RSpace(vec+i) = epsw / MAXVAL(ABS(Blk(1:BlkSiz,i)))
          CSpace(vec+i) = epsw / MAXVAL(ABS(Blk(i,1:BlkSiz)))
        END FORALL
      ELSE
!              {  Rcond <= 1.0D-3  }
!              Condition number of diagonal block exceeds 1.0D3
        nsingBlk = nsingBlk + 1
        
!              Column and Row space remain  0.0D0!
      END IF
      
      vec = vec + BlkSiz
    END DO
    
!        Free workspace for Bd:
    CALL diafree(Bd)
    
!        Free workspace segments:
    DEALLOCATE( temp, STAT=ier)
    IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Deallocation error')
      
      IF (outlev >= 1  .ANd. nsingBlk > 0) THEN
!       Warning, if condition number of diagonal block exceeds 1.0D3
        PRINT '( A, 2X, A, /, 3X,I4, A)' , 'Warning in', rounam, &
        nsingBlk, ' nearly singular block(s) in Block Diagonal matrix!'
      END IF

    ELSE

!     Size of diagonal blocks = 1:

      DO i = 1, S%N
        absval = ABS(S%dia%com(1,i))
        IF (absval > 1.0D-8) THEN
          CSpace(i) = epsw * absval
          RSpace(i) = epsw * absval
        ELSE

!         {  ABS(coAd(i) <= 1.0D-8  }

          nsingBlk = nsingBlk + 1

!         Column and Row space remain  0.0D0!

        END IF
      END DO
      
      IF (outlev >= 1  .ANd. nsingBlk > 0) THEN
        PRINT '( /, A, 2X, A, /, 3X, A, / )' , 'Warning in', rounam,  &
            'Nearly zero diagonal element(s) in Diagonal matrix!'
      END IF
    END IF
    

END SUBROUTINE lumpspace

END MODULE
