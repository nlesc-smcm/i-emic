!     * ================================================================================
!     *
!     * Subroutines for the approximation of the Jacobian Matrix w.r.t.:
!     * - Horizonal mixing of salt and heat according to neutral physics,
!     * - Convective adjustment of salt and heat, 
!     * - Energetically consistant vertical mixing of salt and heat, and 
!     * - Gent-McWilliams stirring of salt and heat,
!     * using subroutines [in mix_sup.f] and method by Coleman et al.
!     *
!     * This file was entirely revised by M. den Toom (2008 => M.dentoom@uu.nl).
!     * The old Cox (1987) implementation proved to be unstable (see Griffies(1998)),
!     * and was therefore replaced by Griffies (1998) alternative. See additional
!     * information in the heading of vmix_fun.
!     *
!     * --------------------------------------------------------------------------------
!     * 
!     * vmix_flag: 
!     * - 0: do nothing
!     * - 1: fixed partition
!     *      ~ evaluated during intialization and remains fixed
!     *      ~ for T,S-equations and T,S-variables, including zero-valued fields
!     * - 2: adaptive partitioning wrt. fields
!     *      ~ partition dependent on T,S-fields, i.e. zero fields are not considered
!     *      ~ checks L_2 norm of the T,S-fields every continuation step and changes
!     *        partition if necessary
!     *
!     * ================================================================================

Call structure DSM and FDJS

!   +----------------------------+                    
!   |                            |   +-----------+    
!   |          mix_imp.f         |   |   usrc.f  |    
!   |                            | --+-----+-----+    
!   |                          --+/        |          
!   |        +---------------+/  |         |          
!   |        | vmix_control  |   |         |          
!   |        +------+--------+   |         |          
!   |               |            +---------+---------+
!   |      +--------+----+         +-------+--+      |
!   |      |  vmix_part  |         | vmix_jac |      |
!   |      +------+---+--+         +----+-----+      |
!   |             |   |                 |            |
!   | +-------------+ |             +---+----+       |
!   | | vmix_el(1,2)| |             |  FDJ   |       |
!   | |  -vmix_dim  | |            /+--------+       |
!   | +-------------+ |           /                  |
!   |                 |          |                   |
!   |             +---+--------+ /                   |
!   |             | DSM        |/                    |
!   |             | -vmix_iptr |                     |
!   |             +------------+                     |
!   |                                                |
!   |                                                |
!   |                                                |
!   +------------------------------------------------+
      
      subroutine vmix_init
      use m_usr
      use m_mix
      implicit none
!include 'usr.com'
!include 'mix.com'

      real time0, time1

      if (vmix.eq.0) then
         vmix_flag = 0
         vmix_diff = 1
         vmix_out  = 1
         vmix_fix  = 1
      else if (vmix.eq.1) then
         vmix_flag = 1
         vmix_diff = 1
         vmix_out  = 1
         vmix_fix  = 1
      else if (vmix.eq.2) then
         vmix_flag = 2
         vmix_diff = 1
         vmix_out  = 1
         vmix_fix  = 0
      else
         vmix_flag = -1
      endif 

      vmix_time=0.0
      call cpu_time(time0)

      if (vmix_out.gt.0) write (99,'(a26)') 'MIX| init...              '
      if (vmix_out.gt.0) write (99,'(a16,i10)') 'MIX|     flag:  ', vmix_flag

      select case (vmix_flag)
      case(0)
         vmix_temp=0
         vmix_salt=0
      case(1)
         vmix_temp=1
         vmix_salt=1
         call vmix_part
      case(2)
         vmix_temp=0
         vmix_salt=0
      end select

      call cpu_time(time1)
      vmix_time=vmix_time+time1-time0
      write (99,'(a26,f10.3)') 'MIX|          ...init done', time1-time0

      end subroutine vmix_init

!     * --------------------------------------------------------------------------------
!     * set mixing parameters
      subroutine vmix_par
      use m_usr
      use m_mix
      implicit none
!include 'usr.com'
!include 'mix.com'

      if (vmix.eq.0) then
         par(MIXP)   =  0.0
         par(P_VC)   =  0.0
         par(ALPC)   =  1.0
         par(ENER)   =  1.0e+2
         par(MKAP)   =  0.0
      endif

      end subroutine vmix_par

!     * --------------------------------------------------------------------------------
      subroutine vmix_control(un)
      use m_usr
      use m_mix
      implicit none

!include 'usr.com'
!include 'mix.com'

      real un(ndim), l2nrm
      integer test_temp, test_salt

      if (vmix_out.gt.0) write (99,'(a26)') 'MIX|   control...         '

      test_temp = 0
      test_salt = 0

! Test for temperature field
      if (l2nrm(un(TT:ndim:nun),n*m*l).gt.1.0e-12) test_temp=1
      if (vmix_out.gt.0)
     +     write (99,'(a16,2i5)') 'MIX|     temp:  ', vmix_temp, test_temp 

! Test for salinity
      if (l2nrm(un(SS:ndim:nun),n*m*l).gt.1.0e-12) test_salt=1
      if (vmix_out.gt.0)
     +     write (99,'(a16,2i5)') 'MIX|     salt:  ', vmix_salt, test_salt 

      if ((vmix_temp.ne.test_temp).or.(vmix_salt.ne.test_salt)) then
         vmix_temp=test_temp 
         vmix_salt=test_salt  
         if ((vmix_temp.ne.0).and.(vmix_temp.ne.0)) call vmix_part
      endif
      vmix_fix=1

      if (vmix_out.gt.0) write (99,'(a26)') 'MIX|       ...control done'

      end subroutine vmix_control
!     * --------------------------------------------------------------------------------
      subroutine vmix_part
      use m_usr
      use m_mix
      implicit none
!include 'usr.com'
!include 'mix.com'

      integer i, info, iwa(6*ndim), liwa
      integer vmix_minrow,vmix_maxrow
      real dnsm

      write (99,'(a26)') 'MIX|   part...            '

      select case(vmix_flag)
      case(1)
         call vmix_el_1(vmix_row, vmix_col, vmix_dim)
      case(2)
         call vmix_el_2(vmix_row, vmix_col, vmix_dim)
      end select

      if (vmix_dim.eq.0) then
         write (*,'(a16,i10)')    'MIX|     idim:  ', vmix_dim
         write (*,*) '   vmix_dim = 0 (NPAIRS), probably no water here'
         write (*,*) '   returning...'
         return
      end if
    
      liwa=6*ndim
      call dsm(ndim,ndim,
     +     vmix_dim,vmix_row,vmix_col,
     +     vmix_ngrp,vmix_maxgrp,vmix_mingrp,
     +     info,
     +     vmix_ipntr,vmix_jpntr,
     +     iwa,liwa)
      dnsm=real(vmix_dim)/(real(ndim)**2)
      write (99,'(a16,i10)')    'MIX|     idim:  ', vmix_dim
      write (99,'(a16,es10.2)') 'MIX|     dnsm:  ', dnsm
      
      if (info.le.0) then
         write (*,*) 'Error in subroutine DSM'
         write (*,*) 'INFO =', info
         write (*,'(a16,i10)')    'MIX|     idim:  ', vmix_dim
         write (*,'(a16,es10.2)') 'MIX|     dnsm:  ', dnsm
         stop
      endif
      vmix_maxrow=0
      vmix_minrow=ndim
      do i=1,ndim
         vmix_maxrow=max(vmix_maxrow,
     +        vmix_ipntr(i+1)-vmix_ipntr(i))
         vmix_minrow=min(vmix_minrow,
     +                  vmix_ipntr(i+1)-vmix_ipntr(i))
      enddo
      write (99,'(a16,i10)') 'MIX|     minrow:', vmix_minrow
      write (99,'(a16,i10)') 'MIX|     maxrow:', vmix_maxrow
      write (99,'(a16,i10)') 'MIX|     mingrp:', vmix_mingrp
      write (99,'(a16,i10)') 'MIX|     maxgrp:', vmix_maxgrp
      write (99,'(a26)')     'MIX|          ...part done'

      end subroutine vmix_part
!     * --------------------------------------------------------------------------------
!     *
!     *=================================================================================
      subroutine vmix_fun(un,mix)

!     *     Computes divergence of the diffusive tracer flux.
!     *     Implementation of neutral physics, Gent-McWilliams stirring and implicit
!     *     vertical mixing. NP and GM are based on Griffies et al. (JPO,1998).
!     *     INPUT
!     *      un:   solution vector
!     *     OUTPUT
!     *      mix: divergence of diffusive tracer flux.
!     *     PARAMETERS
!     *      - MIXP / p6 : isoneutral diffusivity as fraction of PE_H
!     *      - SPL1 / p8 : cut-off value for stability criterion [tprstb]
!     *      - PE_H / p11: horizontal diffusivity
!     *      - PE_V / p12: vertical diffusivity
!     *      - P_VC / p13: implicit vertical diffusivity
!     *      - LAMB / p14: ratio of expansion coefficients
!     *      - NLES / p21: 1.0 for nonlinear equation of state
!     *      - ENER / p24: energy available for consistent mixing
!     *      - ALPC / p25: bifurcation parameter, 0.0 (1.0) is (no) consistent mixing
!     *      - MKAP / p29: GM diffusivity as fraction of PE_H
!     *      - SPL2 / p30: cirtical slope for neutral physics and GM [tprslp]
!     *     TESTING
!     *      - 20/11/08: Established convergence under horizontal grid refinement
!     *                  of order two. Test functions used:
!     *                  {T,S} = \sum_{i=1,3} a_i^{T,S} ( exp(x_i-x_i(0))-1 ).
!     *      - 05/12/08: Experiment (along with earlier tests) confirms that
!     *                  solution is invariant under par(PE_H) when par(MIXP) = 1.0
!     *                  and salinity is constant everywhere, s=0.
!     *
!     *     Author: M. den Toom => m.dentoom@uu.nl
      use m_usr
      use m_mix
      implicit none
!include 'usr.com'
!include 'mix.com'

!     *     Import/export
      real un(ndim), mix(ndim)
!     *     Local
      real u(0:n  ,0:m  ,0:l+1)
      real v(0:n  ,0:m  ,0:l+1)
      real w(0:n+1,0:m+1,0:l  )
      real p(0:n+1,0:m+1,0:l+1)
      real t(0:n+1,0:m+1,0:l+1)
      real s(0:n+1,0:m+1,0:l+1)
      real dtdxe(0:n,0:m+1,0:l+1),dsdxe(0:n,0:m+1,0:l+1)
      real dtdyn(0:n+1,0:m,0:l+1),dsdyn(0:n+1,0:m,0:l+1)
      real dtdzt(0:n+1,0:m+1,0:l),dsdzt(0:n+1,0:m+1,0:l)
      real rho(0:n+1,0:m+1,0:l+1)
      real drhods(0:n+1,0:m+1,0:l+1),drhodt(0:n+1,0:m+1,0:l+1)
      real drhodzt(0:n+1,0:m+1,0:l)
      real Ftxe(0:n,1:m,1:l), Fsxe(0:n,1:m,1:l)
      real Ftyn(1:n,0:m,1:l), Fsyn(1:n,0:m,1:l)
      real Ftzt(1:n,1:m,0:l), Fszt(1:n,1:m,0:l)
      real Ftimp(1:n,1:m,0:l),Fsimp(1:n,1:m,0:l)
      real xes,lambda,piso,pgm,eps,kvc,sp1,sp2
      real dumt,dums,drdh,drdz,slp,tpr
      real, parameter:: epsln = 1.0e-20
      integer i,j,k,ip,jq,kr,row
!     *     Functions
      real tprstb
      integer find_row2

!     *     Extract u,v,w,p,t,s from solution vector un
      call usol(un,u,v,w,p,t,s)      

!     *     Abbreviate the parameter names
      xes    = par(NLES)
      lambda = par(LAMB)
      piso   = par(MIXP)       * par(PE_H)
      pgm    = par(MKAP)       * par(PE_H)
      eps    = (1.0-par(ALPC)) * par(ENER) * par(PE_V) ! as in the old mix_imp.f
      kvc    = par(P_VC)
      sp1    = par(SPL1)
      sp2    = par(SPL2)

!     *     Calculate gradients of temperature and salinity
      call dCdxt(t,dtdxe)       !East face of T-cells
      call dCdxt(s,dsdxe)
      call dCdyt(t,dtdyn)       !North face of T-cells
      call dCdyt(s,dsdyn)
      call dCdzt(t,dtdzt)       !Top face of T-cells
      call dCdzt(s,dsdzt)

!     *     Calculate density and relevant derivatives.
      rho    = lambda*s -  t - xes *
     &     ( alpt1*t +     alpt2*t*t -    alpt3*t*t*t )
      call drhodC(t,drhodt,drhods)
      call dCdzt(rho,drhodzt)

!     *     Calculate fluxes on east, north and top faces ==============================

!     *     Initialize arrays
      Ftxe(:,:,:)   = 0.0
      Fsxe(:,:,:)   = 0.0
      Ftyn(:,:,:)   = 0.0
      Fsyn(:,:,:)   = 0.0
      Ftzt(:,:,:)   = 0.0
      Fszt(:,:,:)   = 0.0
      Ftimp(:,:,:)  = 0.0
      Fsimp(:,:,:)  = 0.0

!     *L0s  start loop over k,j,i
      do k=1,l
         do j=1,m

!     WESTERN BOUNDARY -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
            i=0

!     IFs    NEUTRAL PHYSICS AND GENT-MCWILLIAMS --------------------------------------
            if ( (piso.ne.0.0).or.(pgm.ne.0.0) ) then

!     EAST FACE  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
               dumt   = 0.0
               dums   = 0.0
!     L1s      start loop over 4 triads
               do kr=0,1
                  do ip=0,1
!     Calculate slope and taper
                     drdh   =  ( drhodt(i+ip,j,k)*dtdxe(i,   j,k     ) + 
     &                    drhods(i+ip,j,k)*dsdxe(i,   j,k     ) )
                     drdz   =  ( drhodt(i+ip,j,k)*dtdzt(i+ip,j,k-1+kr) +
     &                    drhods(i+ip,j,k)*dsdzt(i+ip,j,k-1+kr) )
                     call tprslp(drdh,drdz,sp2,slp,tpr)
!     Calculate flux contributions
                     dumt   = dumt + dfzw(k-1+kr) * 
     &                    (  tpr*(piso    )  *     dtdxe(i   ,j,k     ) +
     &                    tpr*(piso-pgm)  * slp*dtdzt(i+ip,j,k-1+kr)  )
                     dums   = dums + dfzw(k-1+kr) * 
     &                    (  tpr*(piso    )  *     dsdxe(i   ,j,k     ) +
     &                    tpr*(piso-pgm)  * slp*dsdzt(i+ip,j,k-1+kr)  )
                  enddo
               enddo
!     L1e      end loop over 4 triads
               Ftxe(i,j,k)  = Ftxe(i,j,k) - dumt/(4*dfzT(k))
               Fsxe(i,j,k)  = Fsxe(i,j,k) - dums/(4*dfzT(k))
            endif
!     IFe    end neutral physics and Gent-McWilliams

!     INTERIOR -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
            do i=1,n

!     IFs    NEUTRAL PHYSICS AND GENT-MCWILLIAMS --------------------------------------
               if ( (piso.ne.0.0).or.(pgm.ne.0.0) ) then

!     EAST FACE  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                  dumt   = 0.0
                  dums   = 0.0
!     L1s      start loop over 4 triads
                  do kr=0,1
                     do ip=0,1
!     Calculate slope and taper
                        drdh   =  ( drhodt(i+ip,j,k)*dtdxe(i,   j,k     ) + 
     &                       drhods(i+ip,j,k)*dsdxe(i,   j,k     ) )
                        drdz   =  ( drhodt(i+ip,j,k)*dtdzt(i+ip,j,k-1+kr) +
     &                       drhods(i+ip,j,k)*dsdzt(i+ip,j,k-1+kr) )
                        call tprslp(drdh,drdz,sp2,slp,tpr)
!     Calculate flux contributions
                        dumt   = dumt + dfzw(k-1+kr) * 
     &                       (  tpr*(piso    )  *     dtdxe(i   ,j,k     ) +
     &                       tpr*(piso-pgm)  * slp*dtdzt(i+ip,j,k-1+kr)  )
                        dums   = dums + dfzw(k-1+kr) * 
     &                       (  tpr*(piso    )  *     dsdxe(i   ,j,k     ) +
     &                       tpr*(piso-pgm)  * slp*dsdzt(i+ip,j,k-1+kr)  )
                     enddo
                  enddo
!     L1e      end loop over 4 triads
                  Ftxe(i,j,k)  = Ftxe(i,j,k) - dumt/(4*dfzT(k))
                  Fsxe(i,j,k)  = Fsxe(i,j,k) - dums/(4*dfzT(k))

!     *         NORTH FACE - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                  dumt   = 0.0
                  dums   = 0.0
!     1s      start loop over 4 triads
                  do kr=0,1
                     do jq=0,1
!     Calculate slope and taper
                        drdh   =  ( drhodt(i,j+jq,k)*dtdyn(i,j   ,k     ) +
     &                       drhods(i,j+jq,k)*dsdyn(i,j   ,k     ) )
                        drdz   =  ( drhodt(i,j+jq,k)*dtdzt(i,j+jq,k-1+kr) +
     &                       drhods(i,j+jq,k)*dsdzt(i,j+jq,k-1+kr) )
                        call tprslp(drdh,drdz,sp2,slp,tpr)
!     Calculate flux contributions
                        dumt   = dumt + dfzw(k-1+kr) * cos(y(j+jq)) *
     &                       (  tpr*(piso    )  *     dtdyn(i,j   ,k     ) +
     &                       tpr*(piso-pgm)  * slp*dtdzt(i,j+jq,k-1+kr)  )
                        dums   = dums + dfzw(k-1+kr) * cos(y(j+jq)) * 
     &                       (  tpr*(piso    )  *     dsdyn(i,   j,k     ) +
     &                       tpr*(piso-pgm)  * slp*dsdzt(i,j+jq,k-1+kr)  )
                     enddo
                  enddo
!     *L1e      end loop over 4 triads
                  Ftyn(i,j,k)  = Ftyn(i,j,k) - dumt/(4*dfzT(k)*cos(yv(j)))
                  Fsyn(i,j,k)  = Fsyn(i,j,k) - dums/(4*dfzT(k)*cos(yv(j)))

!     *         TOP FACE - ZONAL VARIATIONS  - - - -  - - - - - - - - -
                  dumt   = 0.0
                  dums   = 0.0
!     *L1s      start loop over 4 triads
                  do ip=0,1
                     do kr=0,1
!     *         Calculate slope and taper
                        drdh   =  ( drhodt(i,j,k+kr)*dtdxe(i-1+ip,j,k+kr) +
     &                       drhods(i,j,k+kr)*dsdxe(i-1+ip,j,k+kr) )
                        drdz   =  ( drhodt(i,j,k+kr)*dtdzt(i     ,j,k   ) +
     &                       drhods(i,j,k+kr)*dsdzt(i     ,j,k   ) )
                        call tprslp(drdh,drdz,sp2,slp,tpr)
!     *         Calculate flux contributions
                        dumt   = dumt + 
     &                       (  tpr*(piso    ) * slp*slp*dtdzt(i ,j,k   ) +
     &                       tpr*(piso+pgm) * slp*dtdxe(i-1+ip,j,k+kr)  )
                        dums   = dums +
     &                       (  tpr*(piso    ) * slp*slp*dsdzt(i ,j,k   ) +
     &                       tpr*(piso+pgm) * slp*dsdxe(i-1+ip,j,k+kr)  )
                     enddo
                  enddo
!     *L1e      end loop over 4 triads
                  Ftzt(i,j,k)  = Ftzt(i,j,k) - dumt/4
                  Fszt(i,j,k)  = Fszt(i,j,k) - dums/4

!     *         TOP FACE - MERIDIONAL VARIATIONS - - - - - -
                  dumt   = 0.0
                  dums   = 0.0
!     *L1s      start loop over 4 triads
                  do jq=0,1
                     do kr=0,1
!     *         Calculate slope and taper
                        drdh   =  ( drhodt(i,j,k+kr)*dtdyn(i,j-1+jq,k+kr) +
     &                       drhods(i,j,k+kr)*dsdyn(i,j-1+jq,k+kr) )
                        drdz   =  ( drhodt(i,j,k+kr)*dtdzt(i,j     ,k   ) +
     &                       drhods(i,j,k+kr)*dsdzt(i,j     ,k   ) )
                        call tprslp(drdh,drdz,sp2,slp,tpr)
!     *         Calculate flux contributions
                        dumt   = dumt +
     &                       (  tpr*(piso    ) * slp*slp*dtdzt(i ,j,k   ) +
     &                       tpr*(piso+pgm) * slp*dtdyn(i,j-1+jq,k+kr)  )
                        dums   = dums +
     &                       (  tpr*(piso    ) * slp*slp*dsdzt(i ,j,k   ) +
     &                       tpr*(piso+pgm) * slp*dsdyn(i,j-1+jq,k+kr)  )
                     enddo
                  enddo
!     *L1e      end loop over 4 triads
                  Ftzt(i,j,k)  = Ftzt(i,j,k) - dumt/4
                  Fszt(i,j,k)  = Fszt(i,j,k) - dums/4

               endif
!     *IFe    end neutral physics and Gent-McWilliams

!     *IFs    CONSISTENT VERTICAL MIXING -----------------------------------------------
               if (eps.ne.0.0) then
                  Ftzt(i,j,k)  = Ftzt(i,j,k) + tprstb(drhodzt(i,j,k),sp1)*eps *
     &                 dtdzt(i,j,k)/(drhodzt(i,j,k)-epsln)
                  Fszt(i,j,k)  = Fszt(i,j,k) + tprstb(drhodzt(i,j,k),sp1)*eps *
     &                 dsdzt(i,j,k)/(drhodzt(i,j,k)-epsln)
               endif
!     *IFe    end consistent vertical mixing

!     *IFs    IMPLICIT VERTICAL MIXING -------------------------------------------------
               if (kvc.ne.0.0) then
                  Ftimp(i,j,k) = -tprstb(-drhodzt(i,j,k),sp1)*kvc*dtdzt(i,j,k)
                  Fsimp(i,j,k) = -tprstb(-drhodzt(i,j,k),sp1)*kvc*dsdzt(i,j,k)
               endif
!     *IFe    end implicit vertical mixing

            enddo
         enddo
      enddo
!     *L0e  end loop over k,j,i

!     *     Calculate divergence of the fluxes =========================================
!     *L0s  start loop over k,j,i
      do k=1,l
         do j=1,m
            do i=1,n

!     *     For TEMPERATURE ------------------------------------------------------------
               row       = find_row2(i,j,k,TT)
               mix(row)  = 0.0
!     *IFs      
               if ( vmix_temp.eq.1 ) then
!     *     Zonal difference
                  mix(row)  = (Ftxe(i,j,k)           - Ftxe(i-1,j,k)             )/
     &                 (dx*cos(y(j))) + mix(row)
!     *     Meridional difference
                  mix(row)  = (Ftyn(i,j,k)*cos(yv(j))- Ftyn(i,j-1,k)*cos(yv(j-1)))/
     &                 (dy*cos(y(j))) + mix(row)
!     *     Vertical difference
                  mix(row)  = (Ftzt(i,j,k)           - Ftzt(i,j,k-1)             )/
     &                 (dz*dfzT(k))   + mix(row)

                  if (rho_mixing.and.xes.eq.0.0) then
                     mix(row)  = ( (Ftimp(i,j,k) - Ftimp(i,j,k-1))  - 
     &                    (Fsimp(i,j,k) - Fsimp(i,j,k-1))*lambda )/
     &                    ( 2.0*dz*dfzT(k) ) + mix(row)
                  else
                     mix(row)  = (  Ftimp(i,j,k) - Ftimp(i,j,k-1)         )/
     &                    (     dz*dfzT(k) ) + mix(row)
                  endif
               endif
!     *IFe

!     *     For SALINITY ---------------------------------------------------------------
               row       = find_row2(i,j,k,SS)
               mix(row)  = 0.0
!     *IFs
               if ( vmix_salt.eq.1 ) then
!     *     Zonal difference
                  mix(row)  = (Fsxe(i,j,k)           - Fsxe(i-1,j,k)             )/
     &                 (dx*cos(y(j))) + mix(row)
!     *     Meridional difference
                  mix(row)  = (Fsyn(i,j,k)*cos(yv(j))- Fsyn(i,j-1,k)*cos(yv(j-1)))/
     &                 (dy*cos(y(j))) + mix(row)
!     *     Vertical difference
                  mix(row)  = (Fszt(i,j,k)           - Fszt(i,j,k-1)             )/
     &                 (dz*dfzT(k))   + mix(row)
                  if (rho_mixing.and.xes.eq.0.0) then
                     mix(row)  = ( (Fsimp(i,j,k) - Fsimp(i,j,k-1))  - 
     &                    (Ftimp(i,j,k) - Ftimp(i,j,k-1))/lambda )/
     &                    ( 2.0*dz*dfzT(k) ) + mix(row)
                  else
                     mix(row)  = (  Fsimp(i,j,k) - Fsimp(i,j,k-1)         )/
     &                    (     dz*dfzT(k) ) + mix(row)
                  endif
               endif
!     *IFe

            enddo
         enddo
      enddo
!     *L0e  end loop over k,j,i

      end subroutine vmix_fun
!     * --------------------------------------------------------------------------------
      subroutine dCdxt(C,dCdx)

!     *     Calculates zonal derivative of variable that lives on T-grid. Derivative is
!     *     centered on east face of T-cell.
!     *
!     *     Author: M. den Toom => m.dentoom@uu.nl
      use m_usr
      implicit none
!include 'usr.com'

      real    C(0:n+1,0:m+1,0:l+1)
      real dCdx(0:n  ,0:m+1,0:l+1)
      real isoc
      integer i,j,k

      do k=0,l+1
         do j=0,m+1  !--> y might not be defined at 0 and m+1
            do i=0,n
               dCdx(i,j,k) = isoc(i+1,j,k) * isoc(i,j,k) * 
     &              ( C(i+1,j,k) - C(i,j,k) )/( dx * cos(y(j)) )
            enddo
         enddo
      enddo

      end subroutine dCdxt
!     * --------------------------------------------------------------------------------
      subroutine dCdyt(C,dCdy)

!     *     Calculates meridional derivative of variable that lives on T-grid.
!     *     Derivative is centered on north face of T-cell.
!     *
!     *     Author: M. den Toom => m.dentoom@uu.nl
      use m_usr
      implicit none
!include 'usr.com'

      real    C(0:n+1,0:m+1,0:l+1)
      real dCdy(0:n+1,0:m  ,0:l+1)
      real isoc
      integer i,j,k

      do k=0,l+1
         do j=0,m
            do i=0,n+1
               dCdy(i,j,k) = isoc(i,j+1,k) * isoc(i,j,k) * 
     &              ( C(i,j+1,k) - C(i,j,k) )/ dy
            enddo
         enddo
      enddo

      end subroutine dCdyt
!     * --------------------------------------------------------------------------------
      subroutine dCdzt(C,dCdz)

!     *     Calculates vertical derivative of variable that lives on T-grid.
!     *     Derivative is centered on top face of T-cell.
!     *
!     *     Author: M. den Toom => m.dentoom@uu.nl

      use m_usr
      implicit none
!include 'usr.com'

      real    C(0:n+1,0:m+1,0:l+1)
      real dCdz(0:n+1,0:m+1,0:l  )
      real isoc
      integer i,j,k

      do k=0,l
         do j=0,m+1
            do i=0,n+1
               dCdz(i,j,k) = isoc(i,j,k+1) * isoc(i,j,k) * 
     &              ( C(i,j,k+1) - C(i,j,k) )/( dz * dfzW(k) )
            enddo
         enddo
      enddo

      end subroutine dCdzt
!     * --------------------------------------------------------------------------------
      subroutine drhodC(t,drhodt,drhods)

!     *     Calculates expansion coefficients for temperature and salinity.
!     *     Derivative lives on center of T-cell.
!     *
!     *     Author: M. den Toom => m.dentoom@uu.nl

      use m_usr
      implicit none
!include 'usr.com'

      real      t(0:n+1,0:m+1,0:l+1)
      real drhodt(0:n+1,0:m+1,0:l+1),drhods(0:n+1,0:m+1,0:l+1)
      real xes, lambda
      integer i,j,k

!     *     Define ratio of expansion coefficients
      lambda = par(LAMB)
      xes    = par(NLES)

      do k=0,l+1
         do j=0,m+1
            do i=0,n+1
               drhodt(i,j,k) = -1.0 - xes * 
     &              ( alpt1 + 2.0*alpt2*t(i,j,k) - 3.0*alpt3*t(i,j,k)**2 )
               drhods(i,j,k) = lambda
            enddo
         enddo
      enddo

      end subroutine drhodC
!     * --------------------------------------------------------------------------------
      subroutine tprslp(drdh,drdz,spl,slp,tpr)

!     *     Calculates taper tpr based on value of isoneutral slope slp, either
!     *     tap = 1: Gerdes et al. (Clim Dyn, 1991)
!     *     tap = 2: Danabasoglu and McWilliams (J. Clim, 1995)
!     *     tap = 3: De Niet et al. (J. Comp Phys, 2007)
!     *
!     *     Author: M. den Toom => m.dentoom@uu.nl

      use m_usr
      implicit none
!include 'usr.com'

      real drdh,drdz,spl
      real slp,tpr
      real delta,sd,absslp,dum
      real, parameter:: width = 1.0
      real, parameter:: epsln = 1.0e-20

!     *     Calculate slope, cut-off slope, define taper width
      if (drdz.eq.0.0) drdz = epsln
      slp    = -drdh/drdz
      absslp = sqrt( slp**2 ) 
      delta  = (r0dim/hdim)*spl
      sd     = width*delta

!     *     Gerdes et al.
      if     (tap == 1) then 
         if (absslp .gt. delta) then
            tpr = ( delta/absslp )**2
         else
            tpr = 1.0
         endif
!     *     Danabasoglu and McWilliams
      elseif (tap == 2) then
         tpr = 0.5 * (1.0 - tanh( (absslp-delta)/sd  ) )
!     *     De Niet et al
      elseif (tap == 3) then
         if     ( (absslp.lt.delta-sd).and.(drdz.lt.0.0) ) then
            tpr = 1.0
         elseif ( (absslp.ge.delta-sd).and.(absslp.lt.delta)
     &           .and.(drdz.lt.0.0) ) then
            dum = (absslp - (delta-sd))/sd
            tpr = 1.0-3.0*dum**2+2.0*dum**3
         else
            tpr = 0.0
         endif
!     *     Do nothing
      else
         tpr = 1.0
      endif

      end subroutine tprslp
!     * --------------------------------------------------------------------------------
      subroutine vmix_jac(un)
      USE m_mat

!     *     Approximates the jacobian that results from the nonlinear mixing, using the
!     *     method of Coleman et al. (ACM Trans. Math. Software,1984)
!     *
!     *     Adapted by M. den Toom from the old version => m.dentoom@uu.nl

      use m_usr
      use m_mix
      implicit none
!include 'usr.com'
!include 'mix.com'

      real un(ndim),und(ndim),d(ndim)
      real mix(ndim), mixd(ndim)
      real fjac(vmix_dim), eps
      integer i,j,k, numgrp
      integer ix,iy,iz,ie,jx,jy,jz,je,s
      logical col

      eps = 1.0e-08 ! --> adjust?

      select case(vmix_diff)
!     *     Forward differences
      case(1)
         call vmix_fun(un,mix)
         do numgrp=1,vmix_maxgrp
            d = 0.0
            do j=1,ndim
               if (vmix_ngrp(j).eq.numgrp) d(j) = eps
            enddo
            und   = un + d
            call vmix_fun(und, mixd)
            mixd  = mixd - mix
            col   =.false.      ! DO NOT CHANGE
            if (col) then
               call fdjs(ndim,ndim,col,vmix_row,vmix_jpntr,
     +              vmix_ngrp,numgrp,d,mixd,fjac)
            else
               call fdjs(ndim,ndim,col,vmix_col,vmix_ipntr,
     +              vmix_ngrp,numgrp,d,mixd,fjac)
            endif
         enddo
!     *     Central differences      
      case(2)
         do numgrp=1,vmix_maxgrp
            d = 0.0
            do j=1,ndim
               if (vmix_ngrp(j).eq.numgrp) d(j) = 2.0*eps
            enddo
            und   = un + 0.5*d
            call vmix_fun(und, mixd)
            und   = un - 0.5*d
            call vmix_fun(und, mix)
            mixd  = mixd - mix
            col   =.false.      ! DO NOT CHANGE
            if (col) then
               call fdjs(ndim,ndim,col,vmix_row,vmix_jpntr,
     +              vmix_ngrp,numgrp,d,mixd,fjac)
            else
               call fdjs(ndim,ndim,col,vmix_col,vmix_ipntr,
     +              vmix_ngrp,numgrp,d,mixd,fjac)
            endif
         enddo
      end select
! note: row=i, columns are icol(ipntr(i):ipntr(i+1)-1)  

      numgrp = 0
      do i=1,ndim               ! loop over rows, assumes col=.false. !!
         call findex(i,ix,iy,iz,ie)
         do j=vmix_ipntr(i),vmix_ipntr(i+1)-1 ! loop over columns (in approx.)
            call findex(vmix_col(j),jx,jy,jz,je)
            s = 0
            if (jz-iz.eq. -1) s  = s + 9
            if (jz-iz.eq.  1) s  = s + 18 
            if (jx-ix.eq.  0) s  = s + 3
            if (jx-ix.eq.  1) s  = s + 6
            if (jx-ix.eq.1-n) s  = s + 6
            if (jy-iy.eq. -1) s  = s + 1
            if (jy-iy.eq.  0) s  = s + 2
            if (jy-iy.eq.  1) s  = s + 3
            an(s,ie,je,ix,iy,iz) = an(s,ie,je,ix,iy,iz) + fjac(j)
         enddo
      enddo

      end subroutine vmix_jac
!     * --------------------------------------------------------------------------------
      real function isoc(i,j,k)

!     *     Function "is ocean?" (isoc)
!     *
!     *     Author: M. den Toom => m.dentoom@uu.nl
      use m_usr

      implicit none
!include 'usr.com'

      integer i,j,k

      if (landm(i,j,k)==OCEAN .or. landm(i,j,k)==PERIO) then
         isoc = 1.0
      else
         isoc = 0.0
      endif

      end function isoc
!     * --------------------------------------------------------------------------------
      real function tprstb(grad,spl)

!     *     Calculates taper, based on stability (grad = drhodz).
!     *
!     *     Author: M. den Toom => m.dentoom@uu.nl - Adapted from old mix_imp.f

      use m_usr
      implicit none
!include 'usr.com'

      real grad,fac,spl

!     *     Define steepness of taper (note scaling with alphaT to 
!     *     get dimensionless density)
      fac    = alphaT * spl

!     *     Apply taper
      tprstb = max(tanh((-grad*fac)**3),0.0)

      end function tprstb
!     *=================================================================================
!     *
!     * ---------------------------------------------------------------------------- *
      subroutine vmix_el_1(irow,icol,idim) ! for both T and S
      use m_usr
      use m_mix
      implicit none

!include 'usr.com'
!include 'mix.com'

      integer irow(ndim*(nun*np+1)), icol(ndim*(nun*np+1)), idim
      integer p, h, i, j, k, row, tel
      integer find_row2, st(27,4), kind

      irow=0
      icol=0
      idim=0
      p=0

      st( 1,:)=(/ -1, -1, -1, 0 /)
      st( 2,:)=(/ -1,  0, -1, 0 /)
      st( 3,:)=(/ -1,  1, -1, 0 /)
      st( 4,:)=(/  0, -1, -1, 0 /)
      st( 5,:)=(/  0,  0, -1, 1 /)
      st( 6,:)=(/  0,  1, -1, 0 /)
      st( 7,:)=(/  1, -1, -1, 0 /)
      st( 8,:)=(/  1,  0, -1, 0 /)
      st( 9,:)=(/  1,  1, -1, 0 /)
      st(10,:)=(/ -1, -1,  0, 1 /)
      st(11,:)=(/ -1,  0,  0, 1 /)
      st(12,:)=(/ -1,  1,  0, 1 /)
      st(13,:)=(/  0, -1,  0, 1 /)
      st(14,:)=(/  0,  0,  0, 1 /)
      st(15,:)=(/  0,  1,  0, 1 /)
      st(16,:)=(/  1, -1,  0, 1 /)
      st(17,:)=(/  1,  0,  0, 1 /)
      st(18,:)=(/  1,  1,  0, 1 /)
      st(19,:)=(/ -1, -1,  1, 0 /)
      st(20,:)=(/ -1,  0,  1, 0 /)
      st(21,:)=(/ -1,  1,  1, 0 /)
      st(22,:)=(/  0, -1,  1, 0 /)
      st(23,:)=(/  0,  0,  1, 1 /)
      st(24,:)=(/  0,  1,  1, 0 /)
      st(25,:)=(/  1, -1,  1, 0 /)
      st(26,:)=(/  1,  0,  1, 0 /)
      st(27,:)=(/  1,  1,  1, 0 /)

      do k=1,l
         do j=1,m
            do i=1,n
               tel=0
               if (landm(i,j,k).eq.OCEAN) then
                  do h=1,27
                     kind=landm(i+st(h,1),j+st(h,2),k+st(h,3))
                     select case(kind)
                  case(OCEAN)
                     icol(p+tel+1)=find_row2(i+st(h,1),j+st(h,2),k+st(h,3),TT)
                     icol(p+tel+2)=find_row2(i+st(h,1),j+st(h,2),k+st(h,3),SS)
                     tel=tel+2
                  case(PERIO)
                     icol(p+tel+1)=find_row2(i-st(h,1)*(n-1),j+st(h,2),k+st(h,3),TT)
                     icol(p+tel+2)=find_row2(i-st(h,1)*(n-1),j+st(h,2),k+st(h,3),SS)
                     tel=tel+2
                  end select
               enddo
               irow(p+1    :p+tel  )=find_row2(i,j,k,TT)
               irow(p+tel+1:p+2*tel)=find_row2(i,j,k,SS)
               icol(p+tel+1:p+2*tel)=icol(p+1:p+tel)
               p=p+2*tel
            endif
         enddo
      enddo
      enddo
      idim=p
      end subroutine vmix_el_1
!     * ---------------------------------------------------------------------------- *
      subroutine vmix_el_2(irow,icol,idim) ! for T or S
      use m_usr
      use m_mix
      implicit none

!include 'usr.com'
!include 'mix.com'

      integer irow(ndim*(nun*np+1)), icol(ndim*(nun*np+1)), idim
      integer p, h, i, j, k, row, tel
      integer find_row2, st(27,4), kind
      logical flag


      irow=0
      icol=0
      idim=0
      p=0

      st( 1,:)=(/ -1, -1, -1, 0 /)
      st( 2,:)=(/ -1,  0, -1, 0 /)
      st( 3,:)=(/ -1,  1, -1, 0 /)
      st( 4,:)=(/  0, -1, -1, 0 /)
      st( 5,:)=(/  0,  0, -1, 1 /)
      st( 6,:)=(/  0,  1, -1, 0 /)
      st( 7,:)=(/  1, -1, -1, 0 /)
      st( 8,:)=(/  1,  0, -1, 0 /)
      st( 9,:)=(/  1,  1, -1, 0 /)
      st(10,:)=(/ -1, -1,  0, 1 /)
      st(11,:)=(/ -1,  0,  0, 1 /)
      st(12,:)=(/ -1,  1,  0, 1 /)
      st(13,:)=(/  0, -1,  0, 1 /)
      st(14,:)=(/  0,  0,  0, 1 /)
      st(15,:)=(/  0,  1,  0, 1 /)
      st(16,:)=(/  1, -1,  0, 1 /)
      st(17,:)=(/  1,  0,  0, 1 /)
      st(18,:)=(/  1,  1,  0, 1 /)
      st(19,:)=(/ -1, -1,  1, 0 /)
      st(20,:)=(/ -1,  0,  1, 0 /)
      st(21,:)=(/ -1,  1,  1, 0 /)
      st(22,:)=(/  0, -1,  1, 0 /)
      st(23,:)=(/  0,  0,  1, 1 /)
      st(24,:)=(/  0,  1,  1, 0 /)
      st(25,:)=(/  1, -1,  1, 0 /)
      st(26,:)=(/  1,  0,  1, 0 /)
      st(27,:)=(/  1,  1,  1, 0 /)

      do k=1,l
         do j=1,m
            do i=1,n
               row=find_row2(i,j,k,TT)
! Temperature equation ...
               if (vmix_temp.eq.1) then
                  if (landm(i,j,k).eq.OCEAN) then
                     do h=1,27
                        kind=landm(i+st(h,1),j+st(h,2),k+st(h,3))
                        select case(kind)
                     case(OCEAN)
! ... wrt. temperature
                        icol(p+1)=find_row2(i+st(h,1),j+st(h,2),k+st(h,3),TT)
                        irow(p+1)=row
! ... wrt. salinity
                        if (vmix_salt.eq.1) then
                           icol(p+2)=find_row2(i+st(h,1),j+st(h,2),k+st(h,3),SS)
                           irow(p+2)=row
                        endif
                        p=p+1+vmix_salt
                     case(PERIO)
! ... wrt. temperature
                        icol(p+1)=find_row2(i-st(h,1)*(n-1),j+st(h,2),k+st(h,3),TT)
                        irow(p+1)=row
! ... wrt. salinity
                        if (vmix_salt.eq.1) then
                           icol(p+2)=find_row2(i-st(h,1)*(n-1),j+st(h,2),k+st(h,3),SS)
                           irow(p+2)=row
                        endif
                        p=p+1+vmix_salt
                     end select
                  enddo
               endif
            endif
            row = find_row2(i,j,k,SS)
! Salinity equation...
            if (vmix_salt.eq.1) then
               if (landm(i,j,k).eq.OCEAN) then
                  do h=1,27
                     kind=landm(i+st(h,1),j+st(h,2),k+st(h,3))
                     select case(kind)
                  case(OCEAN)
                     if (vmix_temp.eq.1) then
                        icol(p+1)=find_row2(i+st(h,1),j+st(h,2),k+st(h,3),TT)
                        irow(p+1)=row
                     endif
                     icol(p+1+vmix_temp)=
     +                    find_row2(i+st(h,1),j+st(h,2),k+st(h,3),SS)
                     irow(p+1+vmix_temp)=row
                     p=p+1+vmix_temp
                  case(PERIO)
                     if (vmix_temp.eq.1) then
                        icol(p+1)=find_row2(i-st(h,1)*(n-1),j+st(h,2),k+st(h,3),TT)
                        irow(p+1)=row
                     endif
                     icol(p+1+vmix_temp)=
     +                    find_row2(i-st(h,1)*(n-1),j+st(h,2),k+st(h,3),SS)
                     irow(p+1+vmix_temp)=row
                     p=p+1+vmix_temp
                  end select
               enddo
            endif
         endif
      enddo
      enddo
      enddo
      idim=p

      end subroutine vmix_el_2
!     ============================================================================ 
