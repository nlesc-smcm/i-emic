!=======================================================================

PROGRAM nonsym
! Testprogram for the solution of a linear system
!     
!-----------------------------------------------------------------------

#ifndef WITH_UNION

#define csrmatrix  anymatrix

#endif

USE m_build
USE m_matvec
USE m_iniglb
USE m_readpars
USE m_eugen2d

USE m_dump
USE m_wrdstc
USE m_wrdvec
USE m_rdmt
USE m_rdfilenm
USE m_fstrlen
USE m_getunit
USE m_guerrmsg
USE m_rdfilenm
USE m_cmpsolprc

!     Local Parameters:
!     =================

CHARACTER (LEN=*), PARAMETER      :: rounam = 'nonsym'

!     Local Variables:
!     ================

INTEGER 				:: ier
INTEGER 				:: BlkSiz, n
INTEGER 				:: nblock, Mint, probty
DOUBLE PRECISION, DIMENSION(:), POINTER	:: x, b, xRef
TYPE (csrmatrix), POINTER          	:: A
LOGICAL 				:: testx, xGiven
INTEGER 				:: readeq
INTEGER 				:: i
DOUBLE PRECISION 			:: strtch, CVX, CVY
INTEGER 				:: nrsteps
CHARACTER (LEN=80) 			:: basenm
INTEGER 				:: nmlen

#ifdef DEBUG
!     TRACE INFORMATION
PRINT '(A, X, A)' , 'Entry:', rounam
#endif

!     Initialize the global parameter 'OutLev':
CALL iniglb (4)


!     READ THE INPUT DATA and initialise the parameters for the
!     preconditioner:

CALL readpars (readeq, BlkSiz, Mint, probty, strtch, CVX, CVY, testx)

!     Set up the matrix at [A], the right-hand
!     side at [b], and the initial solution in [x]:

IF (readeq == 1) THEN
!        Determine the matrix and an exact solution from the convection-
!        diffusion equation on a 2D-grid:
  CALL eugen2d (probty, Mint, strtch, CVX, CVY, n, A, xRef)
  nrsteps = 1
!      ELSE IF (ReadEq .EQ. 2) THEN
!        Determine the matrix and an exact solution from the convection-
!        diffusion equation on a 3D-grid:
!         CALL eugen3d (ProbTy, Mint, Strtch, CVX, CVY, N, A, xRef, Buf%dbl)
!         NrSteps = 1
ELSE
!        Read input file (base)name containing matrix and RHS:
  CALL rdfilenm ('Enter name data file:', basenm)
  nmlen = fstrlen(basenm)
  
  IF (readeq == 3) THEN
!           Read matrix and right-hand side from ASCII file stored in
!           stencil format. (System from Arthur):
    
!           Read the system and store the matrix in CSR format:
    CALL wrdstc (basenm(1:nmlen), .false., A, b, ier)
    IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Error in wrdstc')
      
      nrsteps = 2
    ELSE IF (readeq == 4) THEN
!     Read matrix and right-hand side from binary files;
!     the matrix is stored in CSR format. (System from Ena)
      
 !     CALL wrdcsr (basenm(1:nmlen) // '.mtr', .true.,  A, ier)
 !     IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Error in wrdcsr')
      CALL rdmt (basenm(1:nmlen) // '.mtr', .true.,  A)
      
      IF (.Not. testx) THEN
        CALL wrdvec (basenm(1:nmlen) // '.rhs', .true., n, b, ier)
        IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Error in wrdvec')
      ENDIF
      nrsteps = 2
    ELSE IF (readeq == 5) THEN
!     Read matrix and right-hand side from ASCII files;
!     the matrix is stored in CSR format. (System from Ena)
       
!      CALL wrdcsr (basenm(1:nmlen) // '.mtr', .false., A, ier)
!      IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Error in wrdcsr')
      CALL rdmt (basenm(1:nmlen) // '.mtr', .false.,  A)
      
      IF (.Not. testx) THEN        
        CALL wrdvec (basenm(1:nmlen) // '.rhs', .false., n, b, ier)
        IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Error in wrdvec')
        xGiven=.FALSE. 
        ALLOCATE ( xRef(1:n), STAT=ier )
        IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
        CALL wrdvec(basenm(1:nmlen) // '.x', .false., n, xRef, ier)
        IF (ier == 0) xGiven=.TRUE.         
      ENDIF   
      nrsteps = 2
    ELSE
      PRINT '(A, I11)' , 'Illegal code to initialise equation:', readeq
      CALL dump(__FILE__,__LINE__,'Illegal code to initialise equation')
    END IF
  END IF
    
  IF ((readeq == 4  .Or. readeq == 5) .And. .Not. testx) THEN
    IF (n /= A%n) THEN
      PRINT '(A, X, A, A, /, 3X, A, 2(/, 3X, A, I11))' ,  &
            'Error detected in', rounam, '!',  &
            'Right hand side does not correspond with matrix!',  &
            'Order of matrix:       ', A%n, 'Length right hand side:', n
      IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Error detected')
    END IF
  END IF
  n=A%n
! Set up the block sizes
  nblock = n/BlkSiz
      
  IF ( MOD(n,nblock) /= 0 ) THEN
    PRINT '(A, I10, 2X, A, I3)' , 'Number of equations', n,  &
          'NOT an integral multiple of block size', BlkSiz
  END IF

  PRINT '(/, 3(A, X, I10,/))' , 'Number of equations:', n,  &
        'Default block size: ', BlkSiz, 'Number of blocks:   ', nblock

  IF ( readeq == 1  .Or. readeq == 2)  THEN
    ALLOCATE( b(1:n), STAT=ier )
    IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
    testx = .true.
  ELSE IF ( testx ) THEN
      ALLOCATE( b(1:n), STAT=ier )
      IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
      ALLOCATE ( xRef(1:n), STAT=ier )
      IF (ier /= 0) CALL dump(__FILE__,__LINE__,'Allocation error')
            
!   Test problem, with exact solution  SIN(1:N)^T:
      FORALL (i = 1:n) xRef(i) = SIN(DBLE(i))

  ELSE
!   Initialise  xRef  with an in range index value:
  END IF
!      testx=.TRUE.
      IF (testx .or. xGiven) THEN
!       Exact solution in 'xRef', determine right hand side in 'b',
!       with   b := A xRef:

        b = 0.0D0
        CALL matvec (n, 1.0D0, csrtoany(A), xRef, b)
      END IF
      IF (xGiven) testx=.TRUE.
      CALL cmpsolprc (BlkSiz, A, b, x, nrsteps, testx, xRef)

!     Normal end:

        END PROGRAM nonsym
