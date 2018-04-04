#include "fdefs.h"

module m_integrals

  use m_par
  use m_usr

  implicit none

contains

  !!------------------------------------------------------------------
  subroutine salt_advection(un, check)

    use, intrinsic :: iso_c_binding

    implicit none
    real(c_double), dimension(ndim), intent(in) :: un
    real(c_double), intent(out) :: check
    integer  i,j,k
    real    u(0:n  ,0:m,0:l+la+1), v(0:n,0:m  ,0:l+la+1)
    real    w(0:n+1,0:m+1,0:l+la  ), p(0:n+1,0:m+1,0:l+la+1)
    real    t(0:n+1,0:m+1,0:l+la+1), s(0:n+1,0:m+1,0:l+la+1)

    call usol(un,u,v,w,p,t,s)

    check= 0.0
    do i = 1, n
       do j = 1, m
          do k = 1, l
             if ( landm(i,j,l) == OCEAN ) then
                check = check + &
                     (u(i,j,k)+u(i,j-1,k))*(s(i+1,j,k)+s(i,j,k))/(4*dx)-  &
                     (u(i-1,j,k)+u(i-1,j-1,k))*(s(i,j,k)+s(i-1,j,k))/(4*dx)+ &
                     (v(i,j,k)+v(i-1,j,k))*(s(i,j+1,k)+s(i,j,k))*cos(yv(j))/(4*dy)- &
                     (v(i,j-1,k)+v(i-1,j-1,k))*(s(i,j,k)+s(i,j-1,k))*cos(yv(j-1))/(4*dy)+ &
                     w(i,j,k)*(s(i,j,k+1)+s(i,j,k))*cos(y(j))/(2*dz*dfzW(k)) - &
                     w(i,j,k-1)*(s(i,j,k)+s(i,j,k-1))*cos(y(j))/(2*dz*dfzW(k-1)) 
             endif
          enddo
       enddo
    enddo

  end subroutine salt_advection

end module m_integrals
