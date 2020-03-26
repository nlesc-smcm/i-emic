#include "fdebug.h"

!   stencil np = 27:
!  +----------++-------++----------+
!  | 12 15 18 || 3 6 9 || 21 24 27 |
!  | 11 14 17 || 2 5 8 || 20 23 26 |
!  | 10 13 16 || 1 4 7 || 19 22 25 |
!  |  below   || center||  above   |
!  +----------++-------++----------+

!    atom(loc,i,j,k)

SUBROUTINE uderiv(type,atom)
  use m_usr
  implicit none
  !     1:  u
  !     2:  uxx
  !     3:  uyy
  !     4:  uzz
  !     IMPORT/EXPORT
  integer type,i,j,k
  real    atom(np,n,m,l)
  real    cosdx2i(0:m),rdy2i,rdz2i,tand2(0:m),cosd2(0:m)
  real    h1,h2
  real    amh,bmh,bmhy
  !
  atom = 0.0
  SELECT CASE(type)
  CASE(1)
     atom(5,:,:,:) = 1.0
  CASE(2)
     ! u_xx
     cosdx2i = (1.0/(cos(yv)*dx))**2
     do j = 1, m - 1
        do i= 1, n
           atom(2,i,j,:) =   amh(yv(j),ih)*cosdx2i(j)
           atom(8,i,j,:) =   amh(yv(j),ih)*cosdx2i(j)
           atom(5,i,j,:) = -(atom(2,i,j,:) + atom(8,i,j,:))
        enddo
     enddo
  CASE(3)
     ! u_yy
     rdy2i = (1.0/dy)**2
     do i=1,n
        do j=1,m-1
           atom(4,i,j,:) = rdy2i * bmh(y(j),ih)*cos(y(j))/cos(yv(j))
           atom(6,i,j,:) = rdy2i * bmh(y(j+1),ih)*cos(y(j+1))/cos(yv(j))
           atom(5,i,j,:) = -(atom(4,i,j,:) + atom(6,i,j,:))
        enddo
     enddo
  CASE(4)
     rdz2i = (1.0/dz)**2
     DO k = 1, l
        h1 = 1./(dfzT(k)*dfzW(k))
        h2 = 1./(dfzT(k)*dfzW(k-1))
        atom(14,:,:,k) = h2*rdz2i
        atom(23,:,:,k) = h1*rdz2i
        atom(5,:,:,k) = -(atom(14,:,:,k) + atom(23,:,:,k))
     ENDDO
  CASE(5)
     tand2 = 1 - tan(yv)*tan(yv)
     DO j = 1, m-1
        atom(5,:,j,:) = bmh(yv(j),ih)*tand2(j)+tan(yv(j))*bmhy(yv(j),ih)
     ENDDO
  CASE(6)
     tand2 = tan(yv)
     cosd2 = cos(yv)
     DO j = 1, m-1
        atom(2,:,j,:)=(bmhy(yv(j),ih)-(amh(yv(j),ih)+bmh(yv(j),ih))*tand2(j))/(dx*cosd2(j))
        atom(8,:,j,:)=-(bmhy(yv(j),ih)-(amh(yv(j),ih)+bmh(yv(j),ih))*tand2(j))/(dx*cosd2(j))
     ENDDO
  END SELECT

end SUBROUTINE uderiv
!****************************************************
SUBROUTINE vderiv(type,atom)
  use m_usr
  implicit none
  !     1:  v
  !     2:  vxx
  !     3:  vyy
  !     4:  vzz
  !     IMPORT/EXPORT
  integer type,i,j,k
  real    atom(np,n,m,l)
  real    cosdx2i(0:m),dy2i,rdz2i,cosd2(0:m),tand2(0:m)
  real    h1,h2
  real    amh,bmh,bmhy
  !
  atom = 0.0
  SELECT CASE(type)
  CASE(1)
     atom(5,:,:,:) = 1.0
  CASE(2)
     ! vxx
     cosdx2i = (1.0/(cos(yv)*dx))**2
     do i=1,n
        do j = 1, m -1
           atom(2,i,j,:) = bmh(yv(j),ih)*cosdx2i(j)
           atom(5,i,j,:) =-2*bmh(yv(j),ih)*cosdx2i(j)
           atom(8,i,j,:) = bmh(yv(j),ih)*cosdx2i(j)
        enddo
     enddo
  CASE(3)
     ! vyy
     dy2i = (1.0/dy)**2
     do i=1,n
        do j=1,m-1
           atom(4,i,j,:) = dy2i* amh(y(j),ih)*cos(y(j))/cos(yv(j))
           atom(6,i,j,:) = dy2i* amh(y(j+1),ih)*cos(y(j+1))/cos(yv(j))
           atom(5,i,j,:) =-(atom(4,i,j,:) + atom(6,i,j,:))
        enddo
     enddo
  CASE(4)
     rdz2i = (1.0/dz)**2
     DO k = 1, l
        h1 = 1./(dfzT(k)*dfzW(k))
        h2 = 1./(dfzT(k)*dfzW(k-1))
        atom(14,:,:,k) = h2*rdz2i
        atom(23,:,:,k) = h1*rdz2i
        atom(5,:,:,k) = -(atom(14,:,:,k) + atom(23,:,:,k))
     ENDDO
  CASE(5)
     DO j = 1, m-1
        atom(5,:,j,:)=bmh(yv(j),ih)-amh(yv(j),ih)*tan(yv(j))*tan(yv(j))+bmhy(yv(j),ih)*tan(yv(j))
     ENDDO
  CASE(6)
     tand2 = tan(yv)
     cosd2 = cos(yv)
     DO j = 1, m-1
        atom(2,:,j,:)=-((amh(yv(j),ih)+bmh(yv(j),ih))*tand2(j)-bmhy(yv(j),ih))/(dx*cosd2(j))
        atom(8,:,j,:)= ((amh(yv(j),ih)+bmh(yv(j),ih))*tand2(j)-bmhy(yv(j),ih))/(dx*cosd2(j))
     ENDDO
  END SELECT

END SUBROUTINE vderiv
!****************************************************
SUBROUTINE pderiv(type,atom)
  use m_usr
  implicit none
  !     1:  ux
  !     2:  vy
  !     3:  wz
  !     IMPORT/EXPORT
  integer type,i,j,k
  real    atom(np,n,m,l)
  real    cos2i(0:m+1),cos2v(0:m),dzi

  !
  atom = 0.0
  SELECT CASE(type)
  CASE(1)
     cos2i = 1.0/(2*cos(y)*dx)
     do j = 1, m
        do i= 1, n
           atom(2,i,j,:) =-cos2i(j)
           atom(4,i,j,:) = cos2i(j)
           atom(1,i,j,:) =-cos2i(j)
           atom(5,i,j,:) = cos2i(j)
        enddo
     enddo
  CASE(2)
     cos2v = cos(yv)
     cos2i = 1./(2*cos(y)*dy)
     do j = 1, m
        do i= 1, n
           atom(4,i,j,:) = - cos2v(j-1) * cos2i(j)
           atom(2,i,j,:) =   cos2v(j) * cos2i(j)
           atom(1,i,j,:) = - cos2v(j-1) * cos2i(j)
           atom(5,i,j,:) =   cos2v(j) * cos2i(j)
        enddo
     enddo
  CASE(3)
     dzi = 1.0/dz
     do k = 1, l
        atom(5,:,:,k) = dzi/dfzT(k)
        atom(14,:,:,k) =-dzi/dfzT(k)
     enddo
  END SELECT

  ! scale continuity equation with constant factor
  !      DO k = 1,l
  !        atom(:,:,:,k) = atom(:,:,:,k) * dfzT(k)
  !      END DO


end SUBROUTINE pderiv
!****************************************************
SUBROUTINE tderiv(type,atom)
  use m_usr
  implicit none
  !     1:  t
  !     2:  t
  !     4:  txx
  !     5:  tyy
  !     6:  tzz
  !     IMPORT/EXPORT
  integer type
  real    atom(np,n,m,l)
  ! LOCAL
  integer i,j,k
  real    cosdx2i(0:m+1),dy2i,dz2i,h1,h2
  !
  atom = 0.0
  SELECT CASE(type)
  CASE(1) ! surface central points
     atom(5,:,:,L) = 1.0
  CASE(2) ! surface central points
     atom(5,:,:,L) = 1.0
  CASE(3)
     cosdx2i = (1.0/(cos(y)*dx))**2
     do i=1,n
        do j=1,m
           do k=1,l
              atom(2,i,j,k) = cosdx2i(j)*(1 - landm(i,j,l))
              atom(5,i,j,k) =-2*cosdx2i(j)*(1 - landm(i,j,l))
              atom(8,i,j,k) = cosdx2i(j)*(1 - landm(i,j,l))
           enddo
        enddo
     enddo
  CASE(4)
     dy2i = (1.0/dy)**2
     do i=1,n
        do k = 1, l
           do j = 1, m
              atom(4,i,j,k) = (dy2i * cos(yv(j-1)) / cos(y(j)))*(1 - landm(i,j,l))
              atom(6,i,j,k) = (dy2i * cos(yv(j)) / cos(y(j)))*(1 - landm(i,j,l))
              atom(5,i,j,k) = -(atom(4,i,j,k) + atom(6,i,j,k))
           enddo
        enddo
     enddo
  CASE(5)
     dz2i = (1.0/dz)**2
     do k=1,l-1
        h1 = 1./(dfzT(k)*dfzW(k))
        h2 = 1./(dfzT(k)*dfzW(k-1))
        do j=1,m
           do i=1,n
              atom(14,i,j,k) = h2*dz2i*(1 - landm(i,j,l))
              atom(23,i,j,k) = h1*dz2i*(1 - landm(i,j,l))
              atom(5,i,j,k) = -(atom(14,i,j,k) + atom(23,i,j,k))
           enddo
        enddo
     enddo
     k = l  !--> boundary condition but not applied boundary.f ???
     h1 = 1./(dfzT(k)*dfzW(k))
     h2 = 1./(dfzT(k)*dfzW(k-1))
     do i=1,n
        do j=1,m
           atom(14,i,j,k) = h2*dz2i*(1 - landm(i,j,l))
           atom(23,i,j,k) = 0.0
           atom(5,i,j,k)  = -(atom(14,i,j,k) + atom(23,i,j,k))
        enddo
     enddo
  CASE(6)
     DO i = 1,n
        DO j = 1,m
           DO k = 1,l
              atom(23,i,j,k) = 1.0 * (1 - landm(i,j,l))
              atom(5,i,j,k) = 1.0 * (1 - landm(i,j,l))
           ENDDO
        ENDDO
     ENDDO
  CASE(7)
     atom(5,:,:,1) = 1.0
  END SELECT

end SUBROUTINE tderiv

!*********************************************************
SUBROUTINE coriolis(type,atom)
  use m_usr
  implicit none
  !     1: fv in the u momentum equation     *   u   u   j+1
  !     2: fu in the v momentum equation     *     v   v    j
  !     1: v_ij + v_i+1j + v_ij-1 + v_i+1j-1 *   u h u h  j
  !     2: u_ij + u_ij+1 + u_i-1j + u_i-1j+1 *     v   v    j-1
  !                                          *     i  i+1
  !     IMPORT/EXPORT                        *  i-1  i
  integer type
  real    atom(np,n,m,l)
  !     LOCAL
  integer i,j
  real      corv(0:m)
  !
  atom = 0.0
  corv = sin(yv) * coriolis_on
  if (type.EQ.1) then
     do i=1,n
        do j=1,m-1
           atom(5,i,j,:) = corv(j)
        enddo
     enddo
  else if (type.EQ.2) then
     do i=1,n
        do j=1,m-1
           atom(5,i,j,:) = corv(j)
        enddo
     enddo
  end if
  !
end SUBROUTINE coriolis

!*********************************************************
SUBROUTINE gradp(type,atom)
  use m_usr
  implicit none
  !     IMPORT/EXPORT
  integer type
  real    atom(np,n,m,l)
  !     LOCAL
  integer i,j,k
  real    cosdxi(0:m),dyi,dzi
  !
  atom = 0.0
  SELECT CASE(type)
  CASE(1)
     cosdxi = 1./(2*cos(yv)*dx)
     do i=1,n
        do j = 1, m -1
           atom(5,i,j,:) =-cosdxi(j)
           atom(6,i,j,:) =-cosdxi(j)
           atom(8,i,j,:) = cosdxi(j)
           atom(9,i,j,:) = cosdxi(j)
        enddo
     enddo
  CASE(2)
     dyi = 1./(2*dy)
     do i=1,n
        do j = 1, m -1
           atom(5,i,j,:) =-dyi
           atom(8,i,j,:) =-dyi
           atom(6,i,j,:) = dyi
           atom(9,i,j,:) = dyi
        enddo
     enddo
  CASE(3)
     dzi = 1./dz
     do k = 1, l
        atom(5,:,:,k) =-dzi/dfzW(k)
        atom(23,:,:,k) = dzi/dfzW(k)
     enddo
  END SELECT

end SUBROUTINE gradp
!*******************************************************
SUBROUTINE masksi(atom, mask)
  use m_usr
  implicit none
  ! sea ice mask atom

  ! IMPORT/EXPORT
  real atom(np,n,m,l)
  real mask(1:n, 1:m)

  atom = 0.0;
  atom(5,:,:,l) = mask;

end SUBROUTINE masksi

!*******************************************************
SUBROUTINE tnlin(type,atom,u,v,w,t)
  use m_usr
  implicit none
  !     nonlinear terms for the t-equation
  !
  !
  !       utx= (u_i-1j+u_ij)*(t_ij+1k-t_ij-1k)
  !       vty= (v_ij+v_ij-1)*(tij+1k-t_ij-1k)
  !       wtz= (w_ijv_ijk-1)*(tijk+1-t_ijk-1)
  !
  !
  !     1:  trT
  !     2:  urTx
  !     3:  Utrx
  !     4:  vrTy
  !     5:  Vtry
  !     6:  wrTz
  !     7:  Wtrz
  !
  ! IMPORT/EXPORT
  integer type
  real    atom(np,n,m,l)
  real    u(0:n  ,0:m,0:l+1),   v(0:n,0:m  ,0:l+1)
  real    w(0:n+1,0:m+1,0:l)
  real    t(0:n+1,0:m+1,0:l+1)
  ! LOCAL
  integer i,j,k,k0,k1
  real    costdxi(0:m+1),tdzi
  ! EXTERNAL
  real lambda, gam, eps

  atom = 0.0
  gam = 1.0e-06
  eps = 1.0
  k0 = 1
  k1 = l
  lambda = par(LAMB)
  !
  SELECT CASE(type)
  CASE(1)                   ! trT
     ! coefficienten voor u met T als basis; hier niet gebruikt.
     atom(5,:,:,:) =1.0
  CASE(2)                   ! urTx
     ! coefficienten voor u met T als basis; hier alleen voor i-1,j (1) en i,j (4)
     costdxi = 1.0/(4*cos(y)*dx)
     DO k = 1, l
        DO j = 1, m
           DO i = 1, n
              atom(2,i,j,k) = -(t(i,j,k)+t(i-1,j,k))*costdxi(j)*(1 - landm(i,j,l))
              atom(4,i,j,k) =  (t(i+1,j,k)+t(i,j,k))*costdxi(j)*(1 - landm(i,j,l))
              atom(1,i,j,k) = -(t(i,j,k)+t(i-1,j,k))*costdxi(j)*(1 - landm(i,j,l))
              atom(5,i,j,k) =  (t(i+1,j,k)+t(i,j,k))*costdxi(j)*(1 - landm(i,j,l))
           ENDDO
        ENDDO
     ENDDO
  CASE(3)                   ! Utrx/(cos y)
     ! coefficienten voor t met U als basis; hier alleen voor i+1,j (7) en i-1,j (1)
     costdxi = 1.0/(4*cos(y)*dx)
     DO k = 1, l
        DO j = 1, m
           DO i = 1, n
              atom(2,i,j,k) = -(u(i-1,j,k)+u(i-1,j-1,k))*costdxi(j)*(1 - landm(i,j,l))
              atom(8,i,j,k) = (u(i,j,k)+u(i,j-1,k))*costdxi(j)*(1 - landm(i,j,l))
              atom(5,i,j,k) = atom(2,i,j,k) + atom(8,i,j,k)
           ENDDO
        ENDDO
     ENDDO
  CASE(4)                   ! vrTy
     ! coefficienten voor v met T als basis; hier alleen voor i,j-1 (3) en i,j (4)
     costdxi = 1.0/(4*cos(y)*dy)
     DO k = 1, l
        DO j = 1, m
           DO i = 1, n
              atom(4,i,j,k) = -costdxi(j)*(t(i,j,k)+t(i,j-1,k))*cos(yv(j-1))*(1 - landm(i,j,l))
              atom(1,i,j,k) = -costdxi(j)*(t(i,j,k)+t(i,j-1,k))*cos(yv(j-1))*(1 - landm(i,j,l))
              atom(5,i,j,k) = costdxi(j)*(t(i,j+1,k)+t(i,j,k))*cos(yv(j))*(1 - landm(i,j,l))
              atom(2,i,j,k) = costdxi(j)*(t(i,j+1,k)+t(i,j,k))*cos(yv(j))*(1 - landm(i,j,l))
           ENDDO
        ENDDO
     ENDDO
  CASE(5)                   ! Vtry
     ! coefficienten voor t met V als basis; hier alleen voor i,j-1 (3) en i,j+1 (5)
     costdxi = 1.0/(4*cos(y)*dy)
     DO k = 1, l
        DO j = 1, m
           DO i = 1, n
              atom(4,i,j,k) = -(v(i,j-1,k)+v(i-1,j-1,k))*costdxi(j)*cos(yv(j-1))*(1 - landm(i,j,l))
              atom(6,i,j,k) = (v(i,j,k)+v(i-1,j,k))*costdxi(j)*cos(yv(j))*(1 - landm(i,j,l))
              atom(5,i,j,k) = atom(4,i,j,k) + atom(6,i,j,k)
           ENDDO
        ENDDO
     ENDDO
     ! coefficienten voor w met T als basis; hier alleen voor i,j,k-1 (3) en i,j,k (4)
  CASE(6)                   ! wrTz
     tdzi = 1.0/(2*dz)
     DO j = 1, m
        DO i = 1, n
           DO k = 1, l-1
              atom(14,i,j,k) = -tdzi*(1 - landm(i,j,l))*(t(i,j,k)+t(i,j,k-1))/dfzT(k)
              atom(5,i,j,k) = tdzi*(1 - landm(i,j,l))*(t(i,j,k+1)+t(i,j,k))/dfzT(k)
           ENDDO
           k = l
           atom(14,i,j,k) = -tdzi*(1 - landm(i,j,l))*(t(i,j,k)+t(i,j,k-1))/dfzT(k)
           atom(5,i,j,k) = 0.0
        ENDDO
     ENDDO
  CASE(7)                   ! Wtrz
     ! coefficienten voor t met W als basis; hier alleen voor i,j,k-1 (8) en i,j,k+1 (9)
     tdzi = 1.0/(2*dz)
     DO k = 1, l
        DO j = 1, m
           DO i = 1, n
              atom(14,i,j,k) = -w(i,j,k-1)*(1 - landm(i,j,l))*tdzi/dfzT(k)
              atom(23,i,j,k) = w(i,j,k)*(1 - landm(i,j,l))*tdzi/dfzT(k)

              atom(5,i,j,k) =  atom(14,i,j,k) + atom(23,i,j,k)
           ENDDO
        ENDDO
     ENDDO

  END SELECT
  !
end SUBROUTINE tnlin
!*****************************************************************
SUBROUTINE wnlin(type,atom,t)
  use m_usr
  implicit none
  !
  !     nonlinear terms for the w-equation
  !
  !     IMPORT/EXPORT
  integer type
  real    atom(np,n,m,l)
  real    t(0:n+1,0:m+1,0:l+1)
  ! LOCAL
  integer i,j,k
  !
  atom = 0.0
  !
  SELECT CASE(type)
  CASE(1)            ! quadratic term jac
     DO k = 1,l-1
        DO j = 1,m
           DO i = 1,n
              atom(23,i,j,k) = (t(i,j,k)+t(i,j,k+1))/2.
              atom(5,i,j,k) = (t(i,j,k)+t(i,j,k+1))/2.
           ENDDO
        ENDDO
     ENDDO
  CASE(2)            ! quadratic term rhs
     DO k = 1,l-1
        DO j = 1,m
           DO i = 1,n
              atom(23,i,j,k) = t(i,j,k+1)/4.
              atom(5,i,j,k) = (t(i,j,k)+2*t(i,j,k+1))/4.
           ENDDO
        ENDDO
     ENDDO
  CASE(3)            ! cubic term jac
     DO k=1,l-1
        DO j = 1,m
           DO i = 1,n
              atom(5,i,j,k) = 0.375*(t(i,j,k)+t(i,j,k+1))**2
              atom(23,i,j,k) = 0.375*(t(i,j,k)+t(i,j,k+1))**2
           ENDDO
        ENDDO
     ENDDO
  CASE(4)            ! cubic term rhs
     DO k=1,l-1
        DO j = 1,m
           DO i = 1,n
              atom(5,i,j,k) = 0.125*(t(i,j,k)*t(i,j,k)+&
                   3*t(i,j,k+1)*t(i,j,k) +&
                   3*t(i,j,k+1)*t(i,j,k+1))
              atom(23,i,j,k) = 0.125*t(i,j,k+1)*t(i,j,k+1)
           ENDDO
        ENDDO
     ENDDO
  END SELECT
  !
END SUBROUTINE wnlin
!*******************************************************
SUBROUTINE unlin(type,atom,u,v,w)
  use m_usr
  implicit none
  !     nonlinear terms for the t-equation
  !
  !     1:  uux
  !     2:  Urux
  !     3:  uvy1
  !     4:  Urvy1
  !     5:  uwz
  !     6:  Urwz
  !     7:  uvy2
  !     8:  Urvy2
  !
  ! IMPORT/EXPORT
  integer type
  real    atom(np,n,m,l)
  real    u(0:n  ,0:m,0:l+1),   v(0:n,0:m  ,0:l+1)
  real    w(0:n+1,0:m+1,0:l  )

  integer i,j,k
  real    costdxi(0:m),tanr(0:m),tdzi(1:l)
  !
  atom = 0.0
  !
  SELECT CASE(type)
  CASE(1)                   ! uux
     costdxi = 1.0/(2*cos(yv)*dx)
     DO j = 1, m
        DO k = 1, l
           DO i = 1, n-1
              atom(8,i,j,k) = u(i+1,j,k)*costdxi(j)
           ENDDO
           DO i = 2,n
              atom(2,i,j,k) = - u(i-1,j,k)*costdxi(j)
           ENDDO
        ENDDO
     ENDDO
  CASE(2)                   ! Urux
     costdxi = 1.0/(2*cos(yv)*dx)
     DO k = 1, l
        DO j = 1, m
           DO i = 1, n-1
              atom(8,i,j,k) = 2*u(i+1,j,k)*costdxi(j)
           ENDDO
           DO i=2,n
              atom(2,i,j,k) = - 2*u(i-1,j,k)*costdxi(j)
           ENDDO
        ENDDO
     ENDDO
  CASE(3)                   ! uvy1
     costdxi = 1.0/(2*cos(yv)*dy)
     DO k = 1, l
        DO i = 1, n
           DO j = 2, m
              atom(4,i,j,k) = -v(i,j-1,k)*cos(yv(j-1))*costdxi(j)
           ENDDO
           DO j=1,m-1
              atom(6,i,j,k) =  v(i,j+1,k)*cos(yv(j+1))*costdxi(j)
           ENDDO
        ENDDO
     ENDDO
  CASE(4)                   ! Urvy1
     costdxi = 1.0/(2*cos(yv)*dy)
     DO k = 1, l
        DO i = 1, n
           DO j = 2, m
              atom(4,i,j,k) =  -u(i,j-1,k)*cos(yv(j-1))*costdxi(j)
           ENDDO
           DO j=1,m-1
              atom(6,i,j,k) =    u(i,j+1,k)*cos(yv(j+1))*costdxi(j)
           ENDDO
        ENDDO
     ENDDO
  CASE(5)                   ! uwz
     tdzi = 1.0/(8*dfzT*dz)
     DO k = 1, l
        DO j = 1, m
           DO i = 1, n
              atom(23,i,j,k) =  (w(i,j,k)+w(i,j+1,k)+w(i+1,j,k)+w(i+1,j+1,k))*tdzi(k)
              atom(14,i,j,k) = -(w(i,j,k-1)+w(i,j+1,k-1)+w(i+1,j,k-1)+w(i+1,j+1,k-1))*tdzi(k)
              atom(5,i,j,k) = atom(14,i,j,k) + atom(23,i,j,k)
           ENDDO
        ENDDO
     ENDDO
  CASE(6)                   ! Urwz
     tdzi = 1.0/(8*dfzT*dz)
     DO j = 1, m
        DO i = 1, n
           DO k = 1, l
              atom(5,i,j,k)  = (u(i,j,k) + u(i,j,k+1))*tdzi(k)
              atom(6,i,j,k)  = (u(i,j,k) + u(i,j,k+1))*tdzi(k)
              atom(8,i,j,k)  = (u(i,j,k) + u(i,j,k+1))*tdzi(k)
              atom(9,i,j,k)  = (u(i,j,k) + u(i,j,k+1))*tdzi(k)
              atom(14,i,j,k) = - (u(i,j,k) + u(i,j,k-1))*tdzi(k)
              atom(15,i,j,k) = - (u(i,j,k) + u(i,j,k-1))*tdzi(k)
              atom(17,i,j,k) = - (u(i,j,k) + u(i,j,k-1))*tdzi(k)
              atom(18,i,j,k) = - (u(i,j,k) + u(i,j,k-1))*tdzi(k)
           ENDDO
        ENDDO
     ENDDO
  CASE(7)                   ! uvy2
     tanr = tan(yv)
     DO k = 1, l
        DO j = 1, m
           DO i = 1, n
              atom(5,i,j,k) = v(i,j,k)*tanr(j)
           ENDDO
        ENDDO
     ENDDO
  CASE(8)                   ! Urvy2
     tanr = tan(yv)
     DO k = 1, l
        DO j = 1, m
           DO i = 1, n
              atom(5,i,j,k) = u(i,j,k)*tanr(j)
           ENDDO
        ENDDO
     ENDDO
  END SELECT
  !
end SUBROUTINE unlin
!*******************************************************
SUBROUTINE vnlin(type,atom,u,v,w)
  use m_usr
  implicit none
  !     nonlinear terms for the t-equation
  !
  !     1:  uvx
  !     2:  uVrx
  !     3:  vvry
  !     4:  Vrvy
  !     5:  vwz
  !     6:  Vrwz
  !     7:  wvrz
  !     8:  Urt2
  !
  ! IMPORT/EXPORT
  integer type
  real    atom(np,n,m,l)
  real    u(0:n  ,0:m,0:l+1),   v(0:n,0:m  ,0:l+1)
  real    w(0:n+1,0:m+1,0:l  )

  integer i,j,k
  real    costdxi(0:m),tanr(0:m),tdzi(1:l)
  !
  atom = 0.0
  !
  SELECT CASE(type)
  CASE(1)                   ! uvx
     costdxi = 1.0/(2*cos(yv)*dx)
     DO k = 1, l
        DO j = 1, m
           DO i = 1, n-1
              atom(8,i,j,k) = u(i+1,j,k)*costdxi(j)
           ENDDO
           DO i=2,n
              atom(2,i,j,k) = - u(i-1,j,k)*costdxi(j)
           ENDDO
        ENDDO
     ENDDO
  CASE(2)                   ! uVrx
     costdxi = 1.0/(2*cos(yv)*dx)
     DO k = 1, l
        DO j = 1, m
           DO i = 1, n-1
              atom(8,i,j,k) = v(i+1,j,k)*costdxi(j)
           ENDDO
           DO i=2,n
              atom(2,i,j,k) = -v(i-1,j,k)*costdxi(j)
           ENDDO
        ENDDO
     ENDDO
  CASE(3)                   ! vvry
     costdxi = 1.0/(2*cos(yv)*dy)
     DO k = 1, l
        DO i = 1, n
           DO j = 1, m-1
              atom(6,i,j,k) =  v(i,j+1,k)*cos(yv(j+1))*costdxi(j)
           ENDDO
           DO j=2,m
              atom(4,i,j,k) = -v(i,j-1,k)*cos(yv(j-1))*costdxi(j)
           ENDDO
        ENDDO
     ENDDO
  CASE(4)                   ! Vrvy
     costdxi = 1.0/(2*cos(yv)*dy)
     DO k = 1, l
        DO i = 1, n
           DO j = 1, m-1
              atom(6,i,j,k) =  2*v(i,j+1,k)*cos(yv(j+1))*costdxi(j)
           ENDDO
           DO j=2,m
              atom(4,i,j,k) =  -2*v(i,j-1,k)*cos(yv(j-1))*costdxi(j)
           ENDDO
        ENDDO
     ENDDO
  CASE(5)                   ! vwz
     tdzi = 1.0/(8*dfzT*dz)
     DO k = 1, l
        DO j = 1, m
           DO i = 1, n
              atom(23,i,j,k) =  (w(i,j,k)+w(i,j+1,k)+w(i+1,j,k)+w(i+1,j+1,k))*tdzi(k)
              atom(14,i,j,k) = -(w(i,j,k-1)+w(i,j+1,k-1)+w(i+1,j,k-1)+w(i+1,j+1,k-1))*tdzi(k)
              atom(5,i,j,k) = atom(14,i,j,k) + atom(23,i,j,k)
           ENDDO
        ENDDO
     ENDDO
  CASE(6)                   ! Vrwz
     tdzi = 1.0/(8*dfzT*dz)
     DO j = 1, m
        DO i = 1, n
           DO k = 1, l
              atom(5,i,j,k) = (v(i,j,k) + v(i,j,k+1))*tdzi(k)
              atom(6,i,j,k) = (v(i,j,k) + v(i,j,k+1))*tdzi(k)
              atom(8,i,j,k) = (v(i,j,k) + v(i,j,k+1))*tdzi(k)
              atom(9,i,j,k) = (v(i,j,k) + v(i,j,k+1))*tdzi(k)
              atom(14,i,j,k) = - (v(i,j,k) + v(i,j,k-1))*tdzi(k)
              atom(15,i,j,k) = - (v(i,j,k) + v(i,j,k-1))*tdzi(k)
              atom(17,i,j,k) = - (v(i,j,k) + v(i,j,k-1))*tdzi(k)
              atom(18,i,j,k) = - (v(i,j,k) + v(i,j,k-1))*tdzi(k)
           ENDDO
        ENDDO
     ENDDO
  CASE(7)                   ! wvrz
     ! coefficienten voor t met W als basis; hier alleen voor i,j,k-1 (8) en i,j,k+1 (9)
     tanr = tan(yv)
     DO k = 1, l
        DO j = 1, m
           DO i = 1, n
              atom(5,i,j,k) = u(i,j,k)*tanr(j)
           ENDDO
        ENDDO
     ENDDO
  CASE(8)                   ! Urt2
     ! coefficienten voor t met W als basis; hier alleen voor i,j,k-1 (8) en i,j,k+1 (9)
     tanr = tan(yv)
     DO k = 1, l
        DO j = 1, m
           DO i = 1, n
              atom(5,i,j,k) = 2*u(i,j,k)*tanr(j)
           ENDDO
        ENDDO
     ENDDO
  END SELECT
  !
end SUBROUTINE vnlin
!***********************************************
real FUNCTION amh(y,ih)
  implicit none
  !     IMPORT/EXPORT
  real     y
  real     ap
  integer ih
  !
  if (ih.eq.0) then
     amh = 1.0
  else
     ap = 10.0
     amh = 1. + ap*exp(-5*y*y)
  endif
  !
END FUNCTION amh
!*************************************************
real FUNCTION bmh(y,ih)
  implicit none
  !     IMPORT/EXPORT
  real     y
  real     bp
  integer ih
  !
  if (ih.eq.0) then
     bmh = 1.0
  else
     bp =10.0
     bmh = 1.0 + bp*exp(-5*y*y)
  endif
  !
END FUNCTION bmh
!***********************************************
real FUNCTION amhy(y,ih)
  implicit none
  !     IMPORT/EXPORT
  real     y
  real     ap
  integer  ih
  !
  if (ih.eq.0) then
     amhy = 0.0
  else
     ap = 10.0
     amhy = -10.*ap*y*exp(-5*y*y)
  endif
  !
END FUNCTION amhy
!*************************************************
real FUNCTION bmhy(y,ih)
  implicit none
  !     IMPORT/EXPORT
  real     y
  real     bp
  integer ih
  !
  if (ih.eq.0) then
     bmhy = 0.0
  else
     bp = 10.0
     bmhy = -10.*bp*y*exp(-5*y*y)
  endif
  !
END FUNCTION bmhy
!**************************************************
!**************************************************
real FUNCTION hs(x1,x2,y1,y2,lab,eps)
  implicit none
  !     IMPORT/EXPORT
  real     x1,x2,y1,y2,lab,eps,h1
  !
  h1 = x2 - x1 - lab*(y2 - y1)
  hs = 0.5*(1 + tanh(h1/eps))
  !
END FUNCTION hs

