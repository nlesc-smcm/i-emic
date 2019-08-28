#include "fdefs.h"

module m_probe

  use m_usr

  implicit none

contains

  !!------------------------------------------------------------------
  subroutine get_atmosphere_t(atmos_temp)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n) :: atmos_temp
    integer :: i,j,pos

    ! Coupling with T (sensible heat flux) would need the atmospheric
    ! temperature tatm.
    if (coupled_T.eq.1) then
       pos = 1
       do j = 1,m
          do i = 1,n
             atmos_temp(pos) = tatm(i,j)
             pos = pos + 1
          end do
       end do
    else
       _INFO_("**WARNING** No coupling, not obtaining any data")
    end if
  end subroutine get_atmosphere_t

  !!------------------------------------------------------------------
  subroutine get_atmosphere_q(atmos_q)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n) :: atmos_q
    integer :: i,j,pos

    ! Coupling with T and S would both need atmospheric humidity qatm
    ! (latent heat flux in T eq., E-P in S eq.).
    if ((coupled_T.eq.1).or.(coupled_S.eq.1)) then
       pos = 1
       do j = 1,m
          do i = 1,n
             atmos_q(pos) = qatm(i,j)
             pos = pos + 1
          end do
       end do
    else
       _INFO_("**WARNING** No coupling, not obtaining any data")
    end if
  end subroutine get_atmosphere_q

  !!------------------------------------------------------------------
  subroutine get_atmosphere_p(atmos_p)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n) :: atmos_p
    integer :: i,j,pos

    ! Coupling with S would need precipitation P (E-P).
    if (coupled_S.eq.1) then
       pos = 1
       do j = 1,m
          do i = 1,n
             atmos_p(pos) = patm(i,j)
             pos = pos + 1
          end do
       end do
    else
       _INFO_("**WARNING** No coupling, not obtaining any data")
    end if
  end subroutine get_atmosphere_p

  !!------------------------------------------------------------------
  subroutine compute_evap(ocean_evap, un)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr
    use m_atm

    implicit none

    ! import/export
    real(c_double), dimension(m*n)  :: ocean_evap
    real(c_double), dimension(ndim) :: un

    ! local
    integer :: i,j,pos, ctr

    ! external
    integer find_row2

    ! reset vector
    ocean_evap = 0

    ctr = 0;
    if ((coupled_T.eq.1).or.(coupled_S.eq.1)) then
       pos = 1
       do j = 1,m
          do i = 1,n
             if (landm(i,j,l).eq.0) then
                ocean_evap(pos) = eo0 + eta * qdim * (((deltat / qdim) * &
                     dqso * un(find_row2(i,j,l,TT)) - &
                     qatm(i,j)))
                ctr = ctr + 1
             endif
             pos = pos + 1

          end do
       end do

       _INFO2_("  thcm probe:  ctr    = ", ctr)
       _INFO2_("  thcm probe:  eta    = ", eta)
       _INFO2_("  thcm probe:  deltat = ", deltat)
       _INFO2_("  thcm probe:  qdim   = ", qdim)
       _INFO2_("  thcm probe:  dqso   = ", dqso)

    else
       _INFO_("**WARNING** No coupling, not obtaining any data")
    end if

  end subroutine compute_evap

  !!------------------------------------------------------------------
  subroutine get_emip(cemip)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n) :: cemip
    integer :: i,j,pos

    pos = 1
    do j = 1,m
       do i = 1,n
          cemip(pos) = emip(i,j) * (1 - landm(i,j,l))
          pos = pos + 1
       end do
    end do
  end subroutine get_emip

  !!------------------------------------------------------------------
  subroutine get_adapted_emip(cemip)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n) :: cemip
    integer :: i,j,pos

    pos = 1
    do j = 1,m
       do i = 1,n
          cemip(pos) = adapted_emip(i,j)
          pos = pos + 1
       end do
    end do
  end subroutine get_adapted_emip

  !!------------------------------------------------------------------
  subroutine get_emip_pert(cspert)

    use, intrinsic :: iso_c_binding
    use m_par  
    use m_usr

    implicit none
    real(c_double), dimension(m*n) :: cspert
    integer :: i,j,pos

    pos = 1
    do j = 1,m
       do i = 1,n
          cspert(pos) = spert(i,j)
          pos = pos + 1
       end do
    end do
  end subroutine get_emip_pert

  !!------------------------------------------------------------------
  subroutine get_salflux(un, salflux, correction, qsoaflux, qsosflux)
    use, intrinsic :: iso_c_binding
    use m_usr
    use m_atm
    use m_ice
    implicit none
    real(c_double), dimension(m*n) :: salflux
    real(c_double), dimension(m*n) :: qsoaflux
    real(c_double), dimension(m*n) :: qsosflux
    real(c_double) :: correction
    integer :: i,j,pos
    real gamma, dedt, pQSnd, QSos, QSoa

    !     IMPORT/EXPORT
    real(c_double),dimension(ndim) ::    un

    !     LOCAL
    real    u(0:n, 0:m,0:l+1), v(0:n,0:m, 0:l+1)
    real    w(0:n+1,0:m+1,0:l), p(0:n+1,0:m+1,0:l+1)
    real    T(0:n+1,0:m+1,0:l+1), S(0:n+1,0:m+1,0:l+1)

    gamma = par(COMB) * par(SALT) 
    dedt  = (deltat / qdim) * dqso

    call usol(un,u,v,w,p,T,S)
    salflux = 0.0

    pQSnd =  par(COMB) * par(SALT) * QSnd
    nus   =  par(COMB) * par(SALT) *  QSnd
    
    pos = 1
    do j = 1,m
       do i = 1,n

          QSos = pQSnd *                        &
               (zeta * (a0 * (s0+S(i,j,l)) - (t0+T(i,j,l))) & ! QTos component
               - ( Qvar * qsa(i,j) + q0 ) )                 & ! QTsa component
               / ( rhodim * Lf )

          ! patm is dimensional, thus we include the background
          ! contributions in this flux as they may not cancel.
          QSoa = pQSnd * ( eo0 + eta * qdim * (     &
               (deltat / qdim) * dqso * T(i,j,l)    &
               - qatm(i,j)) - patm(i,j) )

          qsoaflux(pos) = QSoa / QSnd * (1-msi(i,j))
          qsosflux(pos) = QSos / QSnd * msi(i,j)

          if (gamma.ne.0) then
             if (coupled_S.eq.1) then

                salflux(pos) = (                      &
                     QSoa + msi(i,j) * (QSos - QSoa)  &
                     ) * (1-landm(i,j,l)) / gamma

             else

                salflux(pos) = (1-landm(i,j,l)) *  &
                     (1 - SRES + SRES*par(BIOT)) * emip(i,j) - &
                     SRES * par(BIOT) * S(i,j,l) / gamma

             endif
          endif
          pos = pos + 1
       end do
    end do
    
    call qint(salflux, correction)
    salflux = salflux - correction

    ! The integral correction is stored, for which we take the
    ! continuation parameters into account.
    correction = correction * gamma

  end subroutine get_salflux

  !!------------------------------------------------------------------
  subroutine get_temflux(un, totflux, swflux, shflux, lhflux, siflux, simask)
    use, intrinsic :: iso_c_binding
    use m_par
    use m_usr
    use m_atm
    use m_ice
    
    implicit none
    
    real(c_double), dimension(m*n) :: totflux
    real(c_double), dimension(m*n) :: swflux
    real(c_double), dimension(m*n) :: shflux
    real(c_double), dimension(m*n) :: lhflux
    real(c_double), dimension(m*n) :: siflux
    real(c_double), dimension(m*n) :: simask
    
    integer :: i,j,pos
    real etabi, dedt, dedq, QToa, QSW, QSH, QLH, QTos

    !     IMPORT/EXPORT
    real(c_double),dimension(ndim) :: un

    !     LOCAL
    real  u(0:n  ,0:m,0:l+1),   v(0:n,0:m    ,0:l+1)
    real  w(0:n+1,0:m+1,0:l  ), p(0:n+1,0:m+1,0:l+1)
    real  T(0:n+1,0:m+1,0:l+1), S(0:n+1,0:m+1,0:l+1)

    etabi =  par(COMB) * par(TEMP)

    dedt  =  eta * qdim * (deltat / qdim) * dqso
    ! dedt  = 0.0 !hack
    dedq  = -eta * qdim
    
    call usol(un,u,v,w,p,T,S)

    pos = 1
    do j = 1,m
       do i = 1,n
          if (landm(i,j,l).eq.OCEAN) then
             ! Heat flux forcing from the atmosphere into the ocean.
             ! External and background contributions 
             ! QToa = QSW − QSH − QLH

             QSW = par(COMB) * par(SUNP) *             & ! shortwave heat flux
                  suno(j) * (1-albe0-albed*albe(i,j))

             QSH = Ooa  * (T(i,j,l) - tatm(i,j))         ! sensible heat flux

             QLH = lvsc *                              & ! latent heat flux
                  ( eo0 + dedt*T(i,j,l) + dedq * qatm(i,j) )

             QToa = QSW - QSH - QLH       ! total atmospheric heat flux

             ! Heat flux forcing from sea ice into the ocean. External
             ! and background contributions.
             QTos = QTnd * zeta * (a0 * (s0 + S(i,j,l)) - (t0+T(i,j,l)))

             ! factor out the nondimensionalization and apply seaice mask
             swflux(pos) =  QSW  / QTnd * (1-msi(i,j))
             shflux(pos) = -QSH  / QTnd * (1-msi(i,j))
             lhflux(pos) = -QLH  / QTnd * (1-msi(i,j))
             siflux(pos) =  QTos / QTnd * msi(i,j)
             simask(pos) =  msi(i,j)

             if (coupled_T.eq.0) then
                totflux(pos) = (1 - TRES + TRES*par(BIOT)) * tatm(i,j) - &
                     TRES * par(BIOT) * T(i,j,l) / etabi
             else
                
                totflux(pos) = (1-landm(i,j,l)) *     & !
                     QToa / QTnd + msi(i,j) * (QTos - QToa) / QTnd
                
             endif
          endif
          pos = pos + 1
       end do
    end do

  end subroutine get_temflux


  !!--------------------------------------------------------------------
  ! Obtain shortwave radiation influence field
  subroutine get_suno( sunofield )
    use, intrinsic :: iso_c_binding
    use m_usr
    use m_atm

    real(c_double), dimension(m*n) :: sunofield
    integer :: i,j,pos

    pos = 1;
    do j = 1,m
       do i = 1,n
          sunofield(pos) = suno(j)
          pos = pos + 1
       enddo
    enddo
  end subroutine get_suno

  !!--------------------------------------------------------------------
  ! Return a collection of derivatives of expressions found in lin and
  ! forcing. The sign is taken positive, corresponding to the
  ! implementation in forcing.
  subroutine get_derivatives( un, dftdm, dfsdq, dfsdm, dfsdg)
    use, intrinsic :: iso_c_binding
    use m_usr
    use m_atm
    use m_ice

    implicit none
    real(c_double), dimension(m*n) :: dftdm
    real(c_double), dimension(m*n) :: dfsdq
    real(c_double), dimension(m*n) :: dfsdm
    real(c_double), dimension(m*n) :: dfsdg

    integer :: i,j,pos

    !     IMPORT/EXPORT
    real(c_double),dimension(ndim) :: un

    !     LOCAL
    real    u(0:n  ,0:m,0:l+1), v(0:n,0:m  ,0:l+1)
    real    w(0:n+1,0:m+1,0:l  ), p(0:n+1,0:m+1,0:l+1)
    real    T(0:n+1,0:m+1,0:l+1), S(0:n+1,0:m+1,0:l+1)

    real    QTos, QToa, To, Ta, Ab, So, qa, pa, Ms, qs
    real    QSos, QSoa
    real    pQSnd

    call usol(un,u,v,w,p,T,S)

    dftdm = 0.0 ! temperature equation derivative w.r.t. sea ice mask
    dfsdq = 0.0 ! salinity equation derivative w.r.t. sea ice heat flux
    dfsdm = 0.0 ! salinity equation derivative w.r.t. sea ice mask
    dfsdg = 0.0 ! salinity equation derivative w.r.t. sea ice int corr
    pos = 1
    nus   = par(COMB) * par(SALT) * eta * qdim * QSnd
    pQSnd = par(COMB) * par(SALT) * QSnd
    do j = 1,m
       do i = 1,n
          if (landm(i,j,l).eq.OCEAN) then

             To   = T(i,j,l)     ! sea surface temperature sst
             So   = S(i,j,l)     ! sea surface salinity sss
             Ta   = tatm(i,j)    ! atmosphere temperature
             Ab   = albe(i,j)    ! albedo
             qa   = qatm(i,j)    ! atmosphere humidity
             pa   = patm(i,j)    ! atmosphere precipitation
             Ms   = msi(i,j)     ! sea ice mask
             qs   = qsa(i,j)     ! sea ice heat flux

             if (coupled_T.eq.1)  then

                ! dftdm part ------------------------------

                QTos = QTnd * zeta * (a0 * (So+s0) - (To+t0) ) !

                QToa = par(COMB) * par(SUNP)            &
                     * suno(j) * (1-albe0-albed*Ab)     & ! shortwave heat flux
                     - Ooa * (To - Ta)                  & ! sensible heat flux
                     - lvsc * eta * qdim *              & ! latent heat flux
                     (deltat / qdim * dqso * To - qa)   & 
                     - lvsc * eo0

                dftdm(pos) = QTos - QToa

             endif

             if (coupled_S.eq.1)  then

                ! dfsdq part ------------------------------
                dfsdq(pos) = -pQSnd * Qvar  / (rhodim * Lf) * Ms

                ! dfsdm part ------------------------------
                QSos =   (  &
                     zeta * (a0 * (s0+So) - (t0+To) )    & ! QTos component
                     - ( Qvar * qs + q0 ) )              & ! QTsa component
                     / ( rhodim * Lf )
                
                QSoa = eo0 + eta * qdim * ( &
                     (deltat / qdim) * dqso * To &
                     - qa) - pa 
                                                
                dfsdm(pos) = pQSnd * (QSos - QSoa)
                
                dfsdg(pos) = -1.0;

             endif
          endif
          pos = pos + 1
       enddo
    enddo
    
  end subroutine get_derivatives

end module m_probe
