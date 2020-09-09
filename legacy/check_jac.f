*********************************************************************
      SUBROUTINE check!(un)
      implicit none
      include 'usr.com'
*     LOCAL
      integer i,j,k
      real    un(ndim),fvec(ndim),fvecstar(ndim),del
      real    Achk(ndim,ndim),err(ndim,ndim)
      real    A(ndim,ndim),g05caf
      real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
      real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
      real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)

      open(14,file='errors')
      open(15,file='berrors')
      call init
      call stpnt(un)
      call allocate_mat  ! allocate matrixarrays
      DO i = 1,ndim
       un(i) = 1.0+0.1*g05caf()
      ENDDO
      call usol(un,u,v,w,p,t,s)

      del = 1.0e-06
      DO j=1,ndim
        call rhs(un,fvec)
        un(j) = un(j)+del
        call rhs(un,fvecstar)
        un(j) = un(j)-del
        DO i=1,ndim
          Achk(i,j)= (fvec(i)-fvecstar(i))/del
        ENDDO
      ENDDO
      call matrix(un,0.0)
      call fillsqmat(A)
      err = A-Achk
      DO i=1,ndim
        DO j=1,ndim
          IF (abs(err(i,j)).GT.del) THEN
          write(14,999) i,j,A(i,j),Achk(i,j),err(i,j)
          END IF
          IF (abs(err(i,j)).GT.(1.0e-03)) THEN
          write(15,999) i,j,A(i,j),Achk(i,j),err(i,j)
          END IF
        ENDDO
      ENDDO
*      stop

 999  format(2i4,3e19.6)
      END 
*************************************************************
      SUBROUTINE fillsqmat(A)
      USE m_mat
      implicit none
      include 'usr.com'
!      include 'mat.com'
*      IMPORT/EXPORT
      real     A(ndim,ndim)
      integer  i,v

      A = 0.0
      DO i=1,ndim
       DO v=begA(i),begA(i+1)-1
        IF (abs(coA(v)).GT.1e-10) A(i,jcoA(v)) = coA(v)
       ENDDO
      ENDDO

      END
