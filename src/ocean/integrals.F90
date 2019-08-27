#include "fdefs.h"

module m_integrals

  use m_par
  use m_usr

  implicit none

contains

  !!------------------------------------------------------------------
  ! The interior fluxes should cancel because of the FVM
  ! discretization. At the boundaries, Dirichlet conditions on u,v,w
  ! give zero salt advection fluxes. Hence the total integral should
  ! always be zero + rounding errors.
  subroutine salt_advection(un, check)

    use, intrinsic :: iso_c_binding

    implicit none
    real(c_double), dimension(ndim),   intent(in) :: un
    real(c_double), dimension(n*m*l), intent(out) :: check
    integer  i,j,k,pos
    real    u(0:n  ,0:m,0:l+1),   v(0:n,0:m  ,0:l+1)
    real    w(0:n+1,0:m+1,0:l  ), p(0:n+1,0:m+1,0:l+1)
    real    t(0:n+1,0:m+1,0:l+1), s(0:n+1,0:m+1,0:l+1)

    call usol(un,u,v,w,p,t,s)

    pos = 1
    do k = 1, l
       do j = 1, m   
          do i = 1, n      
             if ( landm(i,j,l) == OCEAN ) then
                check(pos) = &
                     (u(i,j,k)+u(i,j-1,k))*(s(i+1,j,k)+s(i,j,k))/(4*dx)-  &
                     (u(i-1,j,k)+u(i-1,j-1,k))*(s(i,j,k)+s(i-1,j,k))/(4*dx)+ &
                     (v(i,j,k)+v(i-1,j,k))*(s(i,j+1,k)+s(i,j,k))*cos(yv(j))/(4*dy)- &
                     (v(i,j-1,k)+v(i-1,j-1,k))*(s(i,j,k)+s(i,j-1,k))*cos(yv(j-1))/(4*dy)+ &
                     w(i,j,k)*(s(i,j,k+1)+s(i,j,k))*cos(y(j))/(2*dz*dfzW(k)) - &
                     w(i,j,k-1)*(s(i,j,k)+s(i,j,k-1))*cos(y(j))/(2*dz*dfzW(k-1))

             endif
             pos = pos + 1
          enddo
       enddo
    enddo

  end subroutine salt_advection

  !!-------------------------------------------------------------------
  subroutine salt_diffusion(un, check)

    use, intrinsic :: iso_c_binding

    implicit none

    real(c_double), dimension(ndim),   intent(in) :: un
    real(c_double), dimension(n*m*l), intent(out) :: check
    integer  i,j,k,pos
    real    u(0:n  ,0:m,0:l+1),   v(0:n,0:m  ,0:l+1)
    real    w(0:n+1,0:m+1,0:l  ), p(0:n+1,0:m+1,0:l+1)
    real    t(0:n+1,0:m+1,0:l+1), s(0:n+1,0:m+1,0:l+1)
    real    h1, h2, cay, c1, c2

    call usol(un,u,v,w,p,t,s)

    pos = 1
    do k = 1, l
       h1 = 1./(dfzT(k)*dfzW(k))
       h2 = 1./(dfzT(k)*dfzW(k-1))
       do j = 1, m
          cay = cos(y(j))
          c1  = cos(yv(j))
          c2  = cos(yv(j-1))
          do i = 1, n
             if( landm(i,j,k) == OCEAN ) then
                check(pos) =  cos(y(j)) * dfzT(k) * ( &
                     (   s(i+1,j,k)+   s(i-1,j,k)-      2*s(i,j,k))/(dx*dx*cay*cay)+ & 
                     (c1*s(i,j+1,k)+c2*s(i,j-1,k)-(c1+c2)*s(i,j,k))/(dy*dy*cay)    + & 
                     (h1*s(i,j,k+1)+h2*s(i,j,k-1)-(h1+h2)*s(i,j,k))/(dz*dz)        )
             endif
             pos = pos + 1
          enddo
       enddo
    enddo

  end subroutine salt_diffusion

end module m_integrals
