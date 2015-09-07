*     ============================================================================ *
*     
*     Subroutines for the approximation of the Jacobian Matrix w.r.t.:
*     - Horizonal mixing of salt and heat according to neutral physics,
*     - Convective adjustment of salt and heat, 
*     - Energetically consistant vertical mixing of salt and heat, and 
*     - Gent-McWilliams stirring of salt and heat,
*     using subroutines [in mix_sup.f] and method by Coleman et al.
*     
*     ---------------------------------------------------------------------------- *
*     
*     vmix_flag: 
*     - 0: original formulation [obsolete]
*     ~ only convective adjustment and energetically consistant vertical mixing
*     ~ uses analytical expression of the Jacobian obtained with Mathematica
*     - 1: fixed partition
*     ~ evaluated during intialization and remains fixed
*     ~ for T,S-equations and T,S-variables, including zero-valued fields
*     - 2: adaptive partitioning wrt. fields
*     ~ partition dependent on T,S-fields, i.e. zero fields are not considered
*     ~ checks L_2 norm of the T,S-fields every continuation step and changes
*     partition if necessary
*     
*     ---------------------------------------------------------------------------- *
*     
*     Uses the following parameters [see subroutine stpnt in usrc.f]:
*     - MIXP /  6: 
*     - SPL1 /  8: coefficient for vertical mixing and convective adjustment [f_mix]
*     - PE_H / 11:
*     - PE_V / 12:
*     - P_VC / 13: 
*     - LAMB / 14:
*     - ENER / 24: 
*     - ALPC / 25: bifurication parameter for consistant vertical mixing
*     - MKAP / 29: Gent-McWillims stirring (relative to K_H)
*     - SPL2 / 30: coefficient for neutral physics [f_ntr]
*     
*     ============================================================================ *

! Call structure w.r.t DSM and FDJS
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
!   |                                                |
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

      if (vmix_GLB.eq.1) then
         vmix_flag = 2
         vmix_diff = 1
         vmix_out  = 0
         vmix_fix  = 0
      else
         vmix_flag = -1
      endif
      
      vmix_time=0.0
      call cpu_time(time0)

      if (vmix_out.gt.0) write (*,'(a26)') 'MIX| init...              '
      if (vmix_out.gt.0) write (*,'(a16,i10)') 'MIX|     flag:  ', vmix_flag

      select case (vmix_flag)
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
      if (vmix_out.gt.0) write (*,'(a26,f10.3)') 'MIX|          ...init done', time1-time0

      call flush(6)
      end
*     ---------------------------------------------------------------------------- *
      subroutine vmix_par
      use m_usr
      use m_mix
      implicit none
!include 'usr.com'
!include 'mix.com'

      if (vmix_GLB.eq.0) then
         par(MIXP)   =  0.0
         par(P_VC)   =  0.0
         par(ALPC)   =  1.0
         par(ENER)   =  1.0e+2
         par(MKAP)   =  0.0
      endif	
      
      end
*     ---------------------------------------------------------------------------- *
      subroutine vmix_control(un)
      use m_usr
      use m_mix
      implicit none
!include 'usr.com'
!include 'mix.com'
      
      real un(ndim), l2nrm
      integer test_temp, test_salt, i
      if (vmix_out.gt.0) write (*,'(a26)') 'MIX|   control...         '

!     I would like to start with a zero initial condition, so I'm hardcoding this to 1
!      -Erik 
      test_temp = 1  !0
      test_salt = 1  !0
      
! Test for temperature field
      if (l2nrm(un(TT:ndim:nun),n*m*l).gt.1.0e-12) test_temp=1
      if (vmix_out.gt.0)
     +     write (*,'(a16,2i5)') 'MIX|     temp:  ', vmix_temp, test_temp 
      
! Test for salinity
      if (l2nrm(un(SS:ndim:nun),n*m*l).gt.1.0e-12) test_salt=1
      if (vmix_out.gt.0)
     +     write (*,'(a16,2i5)') 'MIX|     salt:  ', vmix_salt, test_salt 
      
      if ((vmix_temp.ne.test_temp).or.(vmix_salt.ne.test_salt)) then
         vmix_temp=test_temp 
         vmix_salt=test_salt  
         if ((vmix_temp.ne.0).and.(vmix_temp.ne.0)) call vmix_part
      endif  
      vmix_fix=1
      
      if (vmix_out.gt.0) write (*,'(a26)') 'MIX|       ...control done'

      end
*     ---------------------------------------------------------------------------- *
      subroutine vmix_part
      use m_usr
      use m_mix
      implicit none
!include 'usr.com'
!include 'mix.com'
      
      integer i, info, iwa(6*ndim), liwa
      integer vmix_minrow,vmix_maxrow
      real dnsm

      if (vmix_out.gt.0) write (*,'(a26)') 'MIX|   part...            '

      select case(vmix_flag)
      case(1)
         call vmix_el_1(vmix_row, vmix_col, vmix_dim)
      case(2)
         call vmix_el_2(vmix_row, vmix_col, vmix_dim)
      end select
      
      liwa=6*ndim

      call dsm(ndim,ndim,
     +     vmix_dim,vmix_row,vmix_col,
     +     vmix_ngrp,vmix_maxgrp,vmix_mingrp,
     +     info,
     +     vmix_ipntr,vmix_jpntr,
     +     iwa,liwa)
      dnsm=real(vmix_dim)/(real(ndim)**2)
      if (vmix_out.eq.0) write (*,'(a16,i10)')    'MIX|     idim:  ', vmix_dim
      if (vmix_out.eq.0) write (*,'(a16,es10.2)') 'MIX|     dnsm:  ', dnsm

      if (info.le.0) then
         write (*,*) 'Error in subroutine DSM'
         write (*,*) 'INFO =', info
         write (*,*) 'NPAIRS (vmix_dim) =', vmix_dim
         if (vmix_dim.eq.0) then
            write(*,*) ' There are no relevant matrix entries,'
            write(*,*) ' our core is probably located on land...'
            return
         else
            stop
         endif
      endif
      vmix_maxrow=0
      vmix_minrow=ndim
      do i=1,ndim
         vmix_maxrow=max(vmix_maxrow,
     +        vmix_ipntr(i+1)-vmix_ipntr(i))
         vmix_minrow=min(vmix_minrow,
     +        vmix_ipntr(i+1)-vmix_ipntr(i))
      enddo
      if (vmix_out.gt.0) then
         write (*,'(a16,i10)') 'MIX|     minrow:', vmix_minrow
         write (*,'(a16,i10)') 'MIX|     maxrow:', vmix_maxrow
         write (*,'(a16,i10)') 'MIX|     mingrp:', vmix_mingrp
         write (*,'(a16,i10)') 'MIX|     maxgrp:', vmix_maxgrp
         write (*,'(a26)')     'MIX|          ...part done'
      endif
      end
*     ---------------------------------------------------------------------------- *
      subroutine vmix_fun(un,mix,mode)
      use m_usr
      use m_mix
      implicit none
!include 'usr.com'
!include 'mix.com'
      real un(ndim), mix(ndim)    ! dezen moet waars geallocate worden
      real u(0:n  ,0:m  ,0:l+la+1)
      real v(0:n  ,0:m  ,0:l+la+1)
      real w(0:n+1,0:m+1,0:l+la  )
      real p(0:n+1,0:m+1,0:l+la+1)
      real t(0:n+1,0:m+1,0:l+la+1)
      real s(0:n+1,0:m+1,0:l+la+1)
      real lambda, kvc, alp, eps0, gmkap, bifm, ph, pv, stkapp
      real T_mmm,T_0mm,T_pmm,T_m0m,T_00m,T_p0m,T_mpm,T_0pm,T_ppm
      real T_mm0,T_0m0,T_pm0,T_m00,T_000,T_p00,T_mp0,T_0p0,T_pp0
      real T_mmp,T_0mp,T_pmp,T_m0p,T_00p,T_p0p,T_mpp,T_0pp,T_ppp
      real S_mmm,S_0mm,S_pmm,S_m0m,S_00m,S_p0m,S_mpm,S_0pm,S_ppm
      real S_mm0,S_0m0,S_pm0,S_m00,S_000,S_p00,S_mp0,S_0p0,S_pp0
      real S_mmp,S_0mp,S_pmp,S_m0p,S_00p,S_p0p,S_mpp,S_0pp,S_ppp
      real dTdx_m00,dTdx_p00,dTdx_0m0,dTdx_0p0,dTdx_00m,dTdx_00p
      real dTdy_m00,dTdy_p00,dTdy_0m0,dTdy_0p0,dTdy_00m,dTdy_00p
      real dTdz_m00,dTdz_p00,dTdz_0m0,dTdz_0p0,dTdz_00m,dTdz_00p
      real dSdx_m00,dSdx_p00,dSdx_0m0,dSdx_0p0,dSdx_00m,dSdx_00p
      real dSdy_m00,dSdy_p00,dSdy_0m0,dSdy_0p0,dSdy_00m,dSdy_00p
      real dSdz_m00,dSdz_p00,dSdz_0m0,dSdz_0p0,dSdz_00m,dSdz_00p
      real dx0, dx1, dy0, dy1, dz2, dz3, dz4, dz6
      real ep, em, kd, f_mix, f_ntr, dum(3)
      real sppk,spmk,smpi,smmi,smpj,smmj,smpk1,smmk1,smpk2,smmk2
      integer i,j,k,row, mode
      integer we, no, es, so, to, bo 
      integer find_row2
      logical flag

      lambda = par(LAMB)
      kvc    = par(P_VC)
      ph     = par(PE_H)
      pv     = par(PE_V)
      alp    = par(ALPC)
      eps0   = par(ENER)
      gmkap  = par(MKAP)
      bifm   = par(MIXP)
      stkapp = par(MKAP)


      call usol(un,u,v,w,p,t,s)      
      mix=0.0
      dum=0.0
      
      do k=1,l
         do j=1,m
            do i=1,n

               if (landm(i,j,k).eq.OCEAN) then
*     MdT  For each OCEAN point define tracer value at that point (000) and the 26
*     neighbouring points; (xxx) are (i,j,k)-pairs, where m=x-1, 0=x, p=x+1.
*     --  landm=PERIO is only relevant for (mxx) and (pxx).
*     --  at all neighbouring points no-flux condition must be implemented. Note that
*     this is problematic for all points with less than two 0's, e.g. T_0mm, because
*     the exact boundary configuration is 'unknown'. I simply substitute the (000)
*     value if such neighbouring point happens to be LAND. (08-04-2008)

! Temperature at stencil points (BC's implemented here!!) ----------------
                  T_mmm=t(i-1,j-1,k-1)*(kd(landm(i-1,j-1,k-1),OCEAN)+kd(landm(i-1,j-1,k-1),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i-1,j-1,k-1),LAND )    
                  T_0mm=t(i  ,j-1,k-1)* kd(landm(i  ,j-1,k-1),OCEAN)
     +                 +t(i  ,j  ,k  )* kd(landm(i  ,j-1,k-1),LAND )
                  T_pmm=t(i+1,j-1,k-1)*(kd(landm(i+1,j-1,k-1),OCEAN)+kd(landm(i+1,j-1,k-1),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i+1,j-1,k-1),LAND )
                  T_m0m=t(i-1,j  ,k-1)*(kd(landm(i-1,j  ,k-1),OCEAN)+kd(landm(i-1,j  ,k-1),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i-1,j  ,k-1),LAND )
                  T_00m=t(i  ,j  ,k-1)* kd(landm(i  ,j  ,k-1),OCEAN)
     +                 +t(i  ,j  ,k  )* kd(landm(i  ,j  ,k-1),LAND )
                  T_p0m=t(i+1,j  ,k-1)*(kd(landm(i+1,j  ,k-1),OCEAN)+kd(landm(i+1,j  ,k-1),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i+1,j  ,k-1),LAND )
                  T_mpm=t(i-1,j+1,k-1)*(kd(landm(i-1,j+1,k-1),OCEAN)+kd(landm(i-1,j+1,k-1),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i-1,j+1,k-1),LAND )
                  T_0pm=t(i  ,j+1,k-1)* kd(landm(i  ,j+1,k-1),OCEAN)
     +                 +t(i  ,j  ,k  )* kd(landm(i  ,j+1,k-1),LAND )
                  T_ppm=t(i+1,j+1,k-1)*(kd(landm(i+1,j+1,k-1),OCEAN)+kd(landm(i+1,j+1,k-1),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i+1,j+1,k-1),LAND )
                  T_mm0=t(i-1,j-1,k  )*(kd(landm(i-1,j-1,k  ),OCEAN)+kd(landm(i-1,j-1,k  ),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i-1,j-1,k  ),LAND )
                  T_0m0=t(i  ,j-1,k  )* kd(landm(i  ,j-1,k  ),OCEAN)
     +                 +t(i  ,j  ,k  )* kd(landm(i  ,j-1,k  ),LAND )
                  T_pm0=t(i+1,j-1,k  )*(kd(landm(i+1,j-1,k  ),OCEAN)+kd(landm(i+1,j-1,k  ),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i+1,j-1,k  ),LAND )
                  T_m00=t(i-1,j  ,k  )*(kd(landm(i-1,j  ,k  ),OCEAN)+kd(landm(i-1,j  ,k  ),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i-1,j  ,k  ),LAND )
                  T_000=t(i  ,j  ,k  )
                  T_p00=t(i+1,j  ,k  )*(kd(landm(i+1,j  ,k  ),OCEAN)+kd(landm(i+1,j  ,k  ),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i+1,j  ,k  ),LAND )
                  T_mp0=t(i-1,j+1,k  )*(kd(landm(i-1,j+1,k  ),OCEAN)+kd(landm(i-1,j+1,k  ),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i-1,j+1,k  ),LAND )
                  T_0p0=t(i  ,j+1,k  )* kd(landm(i  ,j+1,k  ),OCEAN)
     +                 +t(i  ,j  ,k  )* kd(landm(i  ,j+1,k  ),LAND )
                  T_pp0=t(i+1,j+1,k  )*(kd(landm(i+1,j+1,k  ),OCEAN)+kd(landm(i+1,j+1,k  ),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i  ,j+1,k  ),LAND )

                  IF ((landm(i,j,k+1).EQ.ATMOS).OR.(landm(i,j,k+1).EQ.LAND )) THEN
                  T_mmp=t(i-1,j-1,k)
                  T_0mp=t(i  ,j-1,k)
                  T_pmp=t(i+1,j-1,k)
                  T_m0p=t(i-1,j  ,k)
                  T_00p=t(i  ,j  ,k)
                  T_p0p=t(i+1,j  ,k)
                  T_mpp=t(i-1,j+1,k)
                  T_0pp=t(i  ,j+1,k)
                  T_ppp=t(i+1,j+1,k)
               ELSE
                  T_mmp=t(i-1,j-1,k+1)*(kd(landm(i-1,j-1,k+1),OCEAN)+kd(landm(i-1,j-1,k+1),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i-1,j-1,k+1),LAND )
                  T_0mp=t(i  ,j-1,k+1)* kd(landm(i  ,j-1,k+1),OCEAN)
     +                 +t(i  ,j  ,k  )* kd(landm(i  ,j-1,k+1),LAND )
                  T_pmp=t(i+1,j-1,k+1)*(kd(landm(i+1,j-1,k+1),OCEAN)+kd(landm(i+1,j-1,k+1),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i+1,j-1,k+1),LAND )
                  T_m0p=t(i-1,j  ,k+1)*(kd(landm(i-1,j  ,k+1),OCEAN)+kd(landm(i-1,j  ,k+1),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i-1,j  ,k+1),LAND )
                  T_00p=t(i  ,j  ,k+1)* kd(landm(i  ,j  ,k+1),OCEAN)
     +                 +t(i  ,j  ,k  )* kd(landm(i  ,j  ,k+1),LAND )
                  T_p0p=t(i+1,j  ,k+1)*(kd(landm(i+1,j  ,k+1),OCEAN)+kd(landm(i+1,j  ,k+1),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i+1,j  ,k+1),LAND )
                  T_mpp=t(i-1,j+1,k+1)*(kd(landm(i-1,j+1,k+1),OCEAN)+kd(landm(i-1,j+1,k+1),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i-1,j+1,k+1),LAND )
                  T_0pp=t(i  ,j+1,k+1)* kd(landm(i  ,j+1,k+1),OCEAN)
     +                 +t(i  ,j  ,k  )* kd(landm(i  ,j+1,k+1),LAND )
                  T_ppp=t(i+1,j+1,k+1)*(kd(landm(i+1,j+1,k+1),OCEAN)+kd(landm(i+1,j+1,k+1),PERIO))
     +                 +t(i  ,j  ,k  )* kd(landm(i+1,j+1,k+1),LAND )
               ENDIF

! Salinity at stencil points (BC's implemented here!!) -------------------
               S_mmm=s(i-1,j-1,k-1)*(kd(landm(i-1,j-1,k-1),OCEAN)+kd(landm(i-1,j-1,k-1),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i-1,j-1,k-1),LAND )
               S_0mm=s(i  ,j-1,k-1)* kd(landm(i  ,j-1,k-1),OCEAN)
     +              +s(i  ,j  ,k  )* kd(landm(i  ,j-1,k-1),LAND )
               S_pmm=s(i+1,j-1,k-1)*(kd(landm(i+1,j-1,k-1),OCEAN)+kd(landm(i+1,j-1,k-1),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i+1,j-1,k-1),LAND )
               S_m0m=s(i-1,j  ,k-1)*(kd(landm(i-1,j  ,k-1),OCEAN)+kd(landm(i-1,j  ,k-1),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i-1,j  ,k-1),LAND )
               S_00m=s(i  ,j  ,k-1)* kd(landm(i  ,j  ,k-1),OCEAN)
     +              +s(i  ,j  ,k  )* kd(landm(i  ,j  ,k-1),LAND )
               S_p0m=s(i+1,j  ,k-1)*(kd(landm(i+1,j  ,k-1),OCEAN)+kd(landm(i+1,j  ,k-1),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i+1,j  ,k-1),LAND )
               S_mpm=s(i-1,j+1,k-1)*(kd(landm(i-1,j+1,k-1),OCEAN)+kd(landm(i-1,j+1,k-1),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i-1,j+1,k-1),LAND )
               S_0pm=s(i  ,j+1,k-1)* kd(landm(i  ,j+1,k-1),OCEAN)
     +              +s(i  ,j  ,k  )* kd(landm(i  ,j+1,k-1),LAND )
               S_ppm=s(i+1,j+1,k-1)*(kd(landm(i+1,j+1,k-1),OCEAN)+kd(landm(i+1,j+1,k-1),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i+1,j+1,k-1),LAND )
               S_mm0=s(i-1,j-1,k  )*(kd(landm(i-1,j-1,k  ),OCEAN)+kd(landm(i-1,j-1,k  ),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i-1,j-1,k  ),LAND )
               S_0m0=s(i  ,j-1,k  )* kd(landm(i  ,j-1,k  ),OCEAN)
     +              +s(i  ,j  ,k  )* kd(landm(i  ,j-1,k  ),LAND )
               S_pm0=s(i+1,j-1,k  )*(kd(landm(i+1,j-1,k  ),OCEAN)+kd(landm(i+1,j-1,k  ),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i+1,j-1,k  ),LAND )
               S_m00=s(i-1,j  ,k  )*(kd(landm(i-1,j  ,k  ),OCEAN)+kd(landm(i-1,j  ,k  ),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i-1,j  ,k  ),LAND )
               S_000=s(i  ,j  ,k  )
               S_p00=s(i+1,j  ,k  )*(kd(landm(i+1,j  ,k  ),OCEAN)+kd(landm(i+1,j  ,k  ),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i+1,j  ,k  ),LAND )
               S_mp0=s(i-1,j+1,k  )*(kd(landm(i-1,j+1,k  ),OCEAN)+kd(landm(i-1,j+1,k  ),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i-1,j+1,k  ),LAND )
               S_0p0=s(i  ,j+1,k  )* kd(landm(i  ,j+1,k  ),OCEAN)
     +              +s(i  ,j  ,k  )* kd(landm(i  ,j+1,k  ),LAND )
               S_pp0=s(i+1,j+1,k  )*(kd(landm(i+1,j+1,k  ),OCEAN)+kd(landm(i+1,j+1,k  ),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i+1,j+1,k  ),LAND )

               IF ((landm(i,j,k+1).EQ.ATMOS).OR.(landm(i,j,k+1).EQ.LAND )) THEN
               S_mmp=s(i-1,j-1,k)
               S_0mp=s(i  ,j-1,k)
               S_pmp=s(i+1,j-1,k)
               S_m0p=s(i-1,j  ,k)
               S_00p=s(i  ,j  ,k)
               S_p0p=s(i+1,j  ,k)
               S_mpp=s(i-1,j+1,k)
               S_0pp=s(i  ,j+1,k)
               S_ppp=s(i+1,j+1,k)
            ELSE
               S_mmp=s(i-1,j-1,k+1)*(kd(landm(i-1,j-1,k+1),OCEAN)+kd(landm(i-1,j-1,k+1),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i-1,j-1,k+1),LAND )
               S_0mp=s(i  ,j-1,k+1)* kd(landm(i  ,j-1,k+1),OCEAN)
     +              +s(i  ,j  ,k  )* kd(landm(i  ,j-1,k+1),LAND )
               S_pmp=s(i+1,j-1,k+1)*(kd(landm(i+1,j-1,k+1),OCEAN)+kd(landm(i+1,j-1,k+1),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i+1,j-1,k+1),LAND )
               S_m0p=s(i-1,j  ,k+1)*(kd(landm(i-1,j  ,k+1),OCEAN)+kd(landm(i-1,j  ,k+1),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i-1,j  ,k+1),LAND )
               S_00p=s(i  ,j  ,k+1)* kd(landm(i  ,j  ,k+1),OCEAN)
     +              +s(i  ,j  ,k  )* kd(landm(i  ,j  ,k+1),LAND )
               S_p0p=s(i+1,j  ,k+1)*(kd(landm(i+1,j  ,k+1),OCEAN)+kd(landm(i+1,j  ,k+1),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i+1,j  ,k+1),LAND )
               S_mpp=s(i-1,j+1,k+1)*(kd(landm(i-1,j+1,k+1),OCEAN)+kd(landm(i-1,j+1,k+1),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i-1,j+1,k+1),LAND )
               S_0pp=s(i  ,j+1,k+1)* kd(landm(i  ,j+1,k+1),OCEAN)
     +              +s(i  ,j  ,k  )* kd(landm(i  ,j+1,k+1),LAND )
               S_ppp=s(i+1,j+1,k+1)*(kd(landm(i+1,j+1,k+1),OCEAN)+kd(landm(i+1,j+1,k+1),PERIO))
     +              +s(i  ,j  ,k  )* kd(landm(i+1,j+1,k+1),LAND )
            ENDIF

! Derivative of the temperature at the intermediate stencil points -------
            dTdx_m00 = (T_000-T_m00            )/(    dx*cos(y (j  )) ) 
            dTdx_p00 = (T_p00-T_000            )/(    dx*cos(y (j  )) )
            dTdx_0m0 = (T_p00-T_m00+T_pm0-T_mm0)/(4.0*dx*cos(yv(j-1)) )
            dTdx_0p0 = (T_pp0-T_mp0+T_p00-T_m00)/(4.0*dx*cos(yv(j  )) )
            dTdx_00m = (T_p00-T_m00+T_p0m-T_m0m)/(4.0*dx*cos(y (j  )) ) 
            dTdx_00p = (T_p0p-T_m0p+T_p00-T_m00)/(4.0*dx*cos(y (j  )) )

            dTdy_m00 = (T_0p0-T_0m0+T_mp0-T_mm0)/(4.0*dy              )
            dTdy_p00 = (T_pp0-T_pm0+T_0p0-T_0m0)/(4.0*dy              )
            dTdy_0m0 = (T_000-T_0m0            )/(    dy              )
            dTdy_0p0 = (T_0p0-T_000            )/(    dy              )
            dTdy_00m = (T_0p0-T_0m0+T_0pm-T_0mm)/(4.0*dy              )
            dTdy_00p = (T_0pp-T_0mp+T_0p0-T_0m0)/(4.0*dy              )

            dTdz_m00 = (T_00p-T_00m+T_m0p-T_m0m)/(4.0*dz*dfzT(k  )    )
            dTdz_p00 = (T_p0p-T_p0m+T_00p-T_00m)/(4.0*dz*dfzT(k  )    )
            dTdz_0m0 = (T_00p-T_00m+T_0mp-T_0mm)/(4.0*dz*dfzT(k  )    )
            dTdz_0p0 = (T_0pp-T_0pm+T_00p-T_00m)/(4.0*dz*dfzT(k  )    )
            dTdz_00m = (T_000-T_00m            )/(    dz*dfzW(k-1)    )
            dTdz_00p = (T_00p-T_000            )/(    dz*dfzW(k  )    )

! Derivative of the salinity at the intermediate stencil points ----------
            dSdx_m00 = (S_000-S_m00            )/(    dx*cos(y (j  )) )
            dSdx_p00 = (S_p00-S_000            )/(    dx*cos(y (j  )) )
            dSdx_0m0 = (S_p00-S_m00+S_pm0-S_mm0)/(4.0*dx*cos(yv(j-1)) )
            dSdx_0p0 = (S_pp0-S_mp0+S_p00-S_m00)/(4.0*dx*cos(yv(j  )) )
            dSdx_00m = (S_p00-S_m00+S_p0m-S_m0m)/(4.0*dx*cos(y (j  )) )
            dSdx_00p = (S_p0p-S_m0p+S_p00-S_m00)/(4.0*dx*cos(y (j  )) )

            dSdy_m00 = (S_0p0-S_0m0+S_mp0-S_mm0)/(4.0*dy              )
            dSdy_p00 = (S_pp0-S_pm0+S_0p0-S_0m0)/(4.0*dy              )
            dSdy_0m0 = (S_000-S_0m0            )/(    dy              )
            dSdy_0p0 = (S_0p0-S_000            )/(    dy              )
            dSdy_00m = (S_0p0-S_0m0+S_0pm-S_0mm)/(4.0*dy              )
            dSdy_00p = (S_0pp-S_0mp+S_0p0-S_0m0)/(4.0*dy              )

            dSdz_m00 = (S_00p-S_00m+S_m0p-S_m0m)/(4.0*dz*dfzT(k  )    )
            dSdz_p00 = (S_p0p-S_p0m+S_00p-S_00m)/(4.0*dz*dfzT(k  )    )
            dSdz_0m0 = (S_00p-S_00m+S_0mp-S_0mm)/(4.0*dz*dfzT(k  )    )
            dSdz_0p0 = (S_0pp-S_0pm+S_00p-S_00m)/(4.0*dz*dfzT(k  )    )
            dSdz_00m = (S_000-S_00m            )/(    dz*dfzW(k-1)    )
            dSdz_00p = (S_00p-S_000            )/(    dz*dfzW(k)      )

! Functions: 'Heaviside' and energy --------------------------------------
            if (bifm.ne.0.0) then
               smpi =f_ntr((-dTdx_p00+lambda*dSdx_p00),
     +              (-dTdy_p00+lambda*dSdy_p00),
     +              (-dTdz_p00+lambda*dSdz_p00)) ! i+1/2; neutral physics
               smmi =f_ntr((-dTdx_m00+lambda*dSdx_m00),
     +              (-dTdy_m00+lambda*dSdy_m00),
     +              (-dTdz_m00+lambda*dSdz_m00)) ! i-1/2; neutral physics
               smpj =f_ntr((-dTdx_0p0+lambda*dSdx_0p0),
     +              (-dTdy_0p0+lambda*dSdy_0p0),
     +              (-dTdz_0p0+lambda*dSdz_0p0)) ! j+1/2; neutral physics
               smmj =f_ntr((-dTdx_0m0+lambda*dSdx_0m0),
     +              (-dTdy_0m0+lambda*dSdy_0m0),
     +              (-dTdz_0m0+lambda*dSdz_0m0)) ! j-1/2; neutral physics
               smpk2=f_ntr((-dTdx_00p+lambda*dSdx_00p),
     +              (-dTdy_00p+lambda*dSdy_00p),
     +              (-dTdz_00p+lambda*dSdz_00p)) ! k+1/2; neutral physics
               smmk2=f_ntr((-dTdx_00m+lambda*dSdx_00m),
     +              (-dTdy_00m+lambda*dSdy_00m),
     +              (-dTdz_00m+lambda*dSdz_00m)) ! k-1/2; neutral physics
            endif
            if (alp.ne.1.0) then
               smpk1=f_mix(-dTdz_00p+lambda*dSdz_00p) ! k+1/2; stable
               smmk1=f_mix(-dTdz_00m+lambda*dSdz_00m) ! k-1/2; stable
               ep=emix(i,j,k  )
               em=emix(i,j,k-1)
            endif
            sppk =f_mix( dTdz_00p-lambda*dSdz_00p) ! k+1/2; unstable
            spmk =f_mix( dTdz_00m-lambda*dSdz_00m) ! k-1/2; unstable

! dx's, dy's and dz's ----------------------------------------------------
            dx0 = 1/(dx*cos (y(j)) )
            dx1 = 1/(dx*cos (y(j)) )
            dy0 = 1/(dy*cos (y(j)) ) !Note:
            dy1 = 1/(dy*cos (y(j)) ) !  d/dy of divergence also operates on a cos(y)
            dz2 = 1/(dz*dfzT(k   ) )
            dz3 = 1/(dz*dfzT(k   ) )
            dz4 = 1/(dz*dfzT(k   ) )
            dz6 = 1/(dz*dfzT(k   ) )

! Temperature ============================================================ 
            row = find_row2(i,j,k,TT)
            mix(row)=0.0
            
            if (vmix_temp.eq.1) then
               if (bifm.ne.0.0) then
! Temperature Term 0B ----------------------------------------------------
                  if ((1-smpi).gt.0) mix(row)=mix(row)+ph*(1-stkapp)*bifm*
     +                 (1-smpi)*dx0*dTdx_p00
                  if ((1-smmi).gt.0) mix(row)=mix(row)-ph*(1-stkapp)*bifm*
     +                 (1-smmi)*dx0*dTdx_m00
                  if ((1-smpj).gt.0) mix(row)=mix(row)+ph*(1-stkapp)*bifm*
     +                 (1-smpj)*dy0*cos(yv(j))*dTdy_0p0
                  if ((1-smmj).gt.0) mix(row)=mix(row)-ph*(1-stkapp)*bifm*
     +                 (1-smmj)*dy0*cos(yv(j-1))*dTdy_0m0

! Temperature Term 1 -----------------------------------------------------
                  if (smpi.gt.0.0) mix(row)=mix(row)+ph*(1-gmkap)*bifm*
     +                 smpi*dx1*(lambda*dSdx_p00-dTdx_p00)*dTdz_p00/
     +                 (lambda*dSdz_p00-dTdz_p00)
                  if (smmi.gt.0.0) mix(row)=mix(row)-ph*(1-gmkap)*bifm*
     +                 smmi*dx1*(lambda*dSdx_m00-dTdx_m00)*dTdz_m00/
     +                 (lambda*dSdz_m00-dTdz_m00)
                  if (smpj.gt.0.0) mix(row)=mix(row)+ph*(1-gmkap)*bifm*
     +                 smpj*dy1*cos(yv(j))*
     +                 (lambda*dSdy_0p0-dTdy_0p0)*dTdz_0p0/
     +                 (lambda*dSdz_0p0-dTdz_0p0)
                  if (smmj.gt.0.0) mix(row)=mix(row)-ph*(1-gmkap)*bifm*
     +                 smmj*dy1*cos(yv(j-1))*
     +                 (lambda*dSdy_0m0-dTdy_0m0)*dTdz_0m0/
     +                 (lambda*dSdz_0m0-dTdz_0m0)

! Temperature Term 2 -----------------------------------------------------
                  if (smpk2.gt.0.0) mix(row)=mix(row)+ph*(1+gmkap)*bifm*
     +                 smpk2*dz2*((lambda*dSdx_00p*dTdx_00p-dTdx_00p**2)
     +                 + (lambda*dSdy_00p*dTdy_00p-dTdy_00p**2))/
     +                 (lambda*dSdz_00p-dTdz_00p)
                  if (smmk2.gt.0.0) mix(row)=mix(row)-ph*(1+gmkap)*bifm*
     +                 smmk2*dz2*((lambda*dSdx_00m*dTdx_00m-dTdx_00m**2)
     +                 + (lambda*dSdy_00m*dTdy_00m-dTdy_00m**2))/
     +                 (lambda*dSdz_00m-dTdz_00m)

!  Temperature Term 3 ------------------------------------------------------
                  if (smpk2.gt.0.0) mix(row)=mix(row)-ph*bifm*
     +                 smpk2*dz3*(dTdz_00p/((lambda*dSdz_00p-dTdz_00p)**2))*
     +                 ((lambda**2*dSdx_00p**2
     +                 -2.0*lambda*dSdx_00p*dTdx_00p+dTdx_00p**2)
     +                 +(lambda**2*dSdy_00p**2
     +                 -2.0*lambda*dSdy_00p*dTdy_00p+dTdy_00p**2))
                  if (smmk2.gt.0.0) mix(row)=mix(row)+ph*bifm*
     +                 smmk2*dz3*(dTdz_00m/((lambda*dSdz_00m-dTdz_00m)**2))*
     +                 ((lambda**2*dSdx_00m**2
     +                 -2.0*lambda*dSdx_00m*dTdx_00m+dTdx_00m**2)
     +                 +(lambda**2*dSdy_00m**2
     +                 -2.0*lambda*dSdy_00m*dTdy_00m+dTdy_00m**2))
               endif            ! (bifm.gt.0.0)

! Temperature Term 4 -----------------------------------------------------
               if (alp.ne.1.0) then
                  if (smpk1.gt.0.0) mix(row)=mix(row)-(1.0-alp)*eps0*pv*
     +                 dz4*smpk1*dTdz_00p*ep/(dTdz_00p-lambda*dSdz_00p)		
                  if (smmk1.gt.0.0) mix(row)=mix(row)+(1.0-alp)*eps0*pv*
     +                 dz4*smmk1*dTdz_00m*em/(dTdz_00m-lambda*dSdz_00m)
               endif	

! Temperature Term 6 -----------------------------------------------------
               if (kvc.ne.0.0) then
                  if (rho_mixing) then
                     mix(row)= mix(row)-kvc*dz6*
     +                    (sppk*(dTdz_00p-lambda*dSdz_00p)-spmk*(dTdz_00m-lambda*dSdz_00m))/2
                  else
                     mix(row)=mix(row)-kvc*dz6*(sppk*dTdz_00p-spmk*dTdz_00m)
                  endif
               endif
            endif               ! (vmix_temp.eq.1)

! Salinity ===============================================================
            row = find_row2(i,j,k,SS)
            mix(row)=0.0
            
            if (vmix_salt.eq.1) then
               if (bifm.ne.0.0) then          
! Salinity Term 0B ----------------------------------------------------
                  if ((1-smpi).gt.0) mix(row)=mix(row)+ph*(1-stkapp)*bifm*
     +                 (1-smpi)*dx0*dSdx_p00
                  if ((1-smmi).gt.0) mix(row)=mix(row)-ph*(1-stkapp)*bifm*
     +                 (1-smmi)*dx0*dSdx_m00
                  if ((1-smpj).gt.0) mix(row)=mix(row)+ph*(1-stkapp)*bifm*
     +                 (1-smpj)*dy0*cos(yv(j))*dSdy_0p0
                  if ((1-smmj).gt.0) mix(row)=mix(row)-ph*(1-stkapp)*bifm*
     +                 (1-smmj)*dy0*cos(yv(j-1))*dSdy_0m0

! Salinity Term 1 --------------------------------------------------------
                  if (smpi.gt.0.0) mix(row)=mix(row)-ph*(1.0-gmkap)*bifm*
     +                 smpi*dx1*(lambda*dSdx_p00-dTdx_p00)*dSdz_p00/
     +                 (lambda*dSdz_p00-dTdz_p00)
                  if (smmi.gt.0.0) mix(row)=mix(row)+ph*(1.0-gmkap)*bifm*
     +                 smmi*dx1*(lambda*dSdx_m00-dTdx_m00)*dSdz_m00/
     +                 (lambda*dSdz_m00-dTdz_m00)
                  if (smpj.gt.0.0) mix(row)=mix(row)-ph*(1.0-gmkap)*bifm*
     +                 smpj*dy1*cos(yv(j))*
     +                 (lambda*dSdy_0p0-dTdy_0p0)*dSdz_0p0/
     +                 (lambda*dSdz_0p0-dTdz_0p0)
                  if (smmj.gt.0.0) mix(row)=mix(row)+ph*(1.0-gmkap)*bifm*
     +                 smmj*dy1*cos(yv(j-1))*
     +                 (lambda*dSdy_0m0-dTdy_0m0)*dSdz_0m0/
     +                 (lambda*dSdz_0m0-dTdz_0m0)

! Salinity Term 2 --------------------------------------------------------
                  if (smpk2.gt.0.0) mix(row)=mix(row)-ph*(1.0+gmkap)*bifm*
     +                 smpk2*dz2*((lambda*dSdx_00p**2-dSdx_00p*dTdx_00p**2)
     +                 + (lambda*dSdy_00p**2-dSdy_00p*dTdy_00p**2))/
     +                 (lambda*dSdz_00p-dTdz_00p)
                  if (smmk2.gt.0.0) mix(row)=mix(row)+ph*(1.0+gmkap)*bifm*
     +                 smmk2*dz2*((lambda*dSdx_00m**2-dSdx_00m*dTdx_00m**2)
     +                 + (lambda*dSdy_00m**2-dSdy_00m*dTdy_00m**2))/
     +                 (lambda*dSdz_00m-dTdz_00m)

! Salinity Term 3 --------------------------------------------------------
                  if (smpk2.gt.0.0) mix(row)=mix(row)-ph*bifm*
     +                 smpk2**2*dz3*(dSdz_00p/((lambda*dSdz_00p-dTdz_00p)**2))*
     +                 ((lambda**2*dSdx_00p**2
     +                 -2.0*lambda*dSdx_00p*dTdx_00p+dTdx_00p**2)
     +                 +(lambda**2*dSdy_00p**2
     +                 -2.0*lambda*dSdy_00p*dTdy_00p+dTdy_00p**2))
                  if (smmk2.gt.0.0) mix(row)=mix(row)+ph*bifm*
     +                 smmk2**2*dz3*(dSdz_00m/((lambda*dSdz_00m-dTdz_00m)**2))*
     +                 ((lambda**2*dSdx_00m**2
     +                 -2.0*lambda*dSdx_00m*dTdx_00m+dTdx_00m**2)
     +                 +(lambda**2*dSdy_00m**2
     +                 -2.0*lambda*dSdy_00m*dTdy_00m+dTdy_00m**2))
               endif            ! (bifm.gt.0.0)
               
! Salinity Term 4 --------------------------------------------------------
               if (alp.ne.1.0) then
                  if (smpk1.gt.0.0) mix(row)=mix(row)-(1.0-alp)*eps0*pv*
     +                 dz4*smpk1*dSdz_00p*ep/(dTdz_00p-lambda*dSdz_00p)
                  if (smmk1.gt.0.0) mix(row)=mix(row)+(1.0-alp)*eps0*pv*
     +                 dz4*smmk1*dSdz_00m*em/(dTdz_00m-lambda*dSdz_00m)
               endif

! Salinity Term 6 --------------------------------------------------------
               if (kvc.ne.0.0) then
                  if (rho_mixing) then
                     mix(row)=mix(row)-kvc*dz6*
     +                    (sppk*(dSdz_00p-dTdz_00p/lambda)-spmk*(dSdz_00m-dTdz_00m/lambda))/2
                  else
                     mix(row)=mix(row)-kvc*dz6*(sppk*dSdz_00p-spmk*dSdz_00m)
                  endif
               endif
            endif               ! (vmix_salt.eq.1)
            
            if (mode.eq.0) then
               flag=.false.
               if (smpi .gt.0.0) flag=.true.
               if (smmi .gt.0.0) flag=.true.
               if (smpj .gt.0.0) flag=.true.
               if (smmj .gt.0.0) flag=.true.
               if (smpk2.gt.0.0) flag=.true.
               if (smmk2.gt.0.0) flag=.true.
               if (flag) dum(1)=dum(1)+1.0
               flag=.false.
               if (smpk1.gt.0.0) flag=.true.
               if (smmk1.gt.0.0) flag=.true.
               if (flag) dum(2)=dum(2)+1.0
               flag=.false.
               if (sppk .gt.0.0) flag=.true.
               if (spmk .gt.0.0) flag=.true.
               if (flag) dum(3)=dum(3)+1.0
            endif
         endif	    
      enddo
      enddo
      enddo
      
      if(mode.eq.0) then
         dum(1)=100.0*dum(1)/real(n*m*l) ! ntrl. phy.
         dum(2)=100.0*dum(2)/real(n*m*l) ! cons. mix.
         dum(3)=100.0*dum(3)/real(n*m*l) ! conv. adj.
         if (vmix_out.gt.0) write(*,'(a16,3f7.2)') 'MIX|     fun:   ', dum
      endif
      
      end
*     ---------------------------------------------------------------------------- *
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
      end
*     ---------------------------------------------------------------------------- *
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

      end
*     ---------------------------------------------------------------------------- *
      subroutine vmix_jac(un)
      
      use m_mat
      use m_usr
      use m_mix
      implicit none
!include 'usr.com'
!include 'mix.com'
      
      integer i,j,k, numgrp
      integer ix,iy,iz,ie,jx,jy,jz,je,s
      real un(ndim), mix(ndim), d(ndim), und(ndim), mixd(ndim)
      real fjac(vmix_dim), eps
      logical col


      eps=1.0e-08

      select case(vmix_diff)
      case(1)                   ! forward diffences
      
         call vmix_fun(un,mix,1)
         do numgrp=1,vmix_maxgrp
            d=0.0
            do j=1,ndim
               if (vmix_ngrp(j).eq.numgrp) d(j) = eps
            enddo
            und=un+d
            call vmix_fun(und, mixd,1)
            mixd=mixd-mix
            col=.false.         ! DO NOT CHANGE
            if (col) then
               call fdjs(ndim,ndim,col,vmix_row,vmix_jpntr,
     +              vmix_ngrp,numgrp,d,mixd,fjac)
            else
               call fdjs(ndim,ndim,col,vmix_col,vmix_ipntr,
     +              vmix_ngrp,numgrp,d,mixd,fjac)
            endif
         enddo
         
      case(2)                   ! central differences
         do numgrp=1,vmix_maxgrp
            d=0.0
            do j=1,ndim
               if (vmix_ngrp(j).eq.numgrp) d(j) = 2.0*eps
            enddo
            und=un+0.5*d
            call vmix_fun(und, mixd,1)
            und=un-0.5*d
            call vmix_fun(und, mix,1)
            mixd=mixd-mix
            col=.false.         ! DO NOT CHANGE
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
      
      numgrp=0
      do i=1,ndim               ! loop over rows, assumes col=.false. !!
         call findex(i,ix,iy,iz,ie)
         do j=vmix_ipntr(i),vmix_ipntr(i+1)-1 ! loop over columns (in approx.)
            call findex(vmix_col(j),jx,jy,jz,je)
            s=0
            if (jz-iz.eq. -1) s=s+9
            if (jz-iz.eq.  1) s=s+18	  
            if (jx-ix.eq.  0) s=s+3
            if (jx-ix.eq.  1) s=s+6
            if (jx-ix.eq.1-n) s=s+6
            if (jy-iy.eq. -1) s=s+1
            if (jy-iy.eq.  0) s=s+2
            if (jy-iy.eq.  1) s=s+3
            al(ix,iy,iz,s,ie,je)=al(ix,iy,iz,s,ie,je)+fjac(j)
         enddo
      enddo

      end
*     ---------------------------------------------------------------------------- *
      real function f_mix(val)
      use m_usr
      implicit none
      !include 'usr.com'
      
      real val, fac
      
      fac=par(SPL1)
!     if (val.gt.0.0) then 
!     f_mix = 0.0      
!     else
!     f_mix=max(tanh((val*fac)**2),0.0)
!     endif
      f_mix=max(tanh(-(val*fac)**3),0.0)

      end
*     ---------------------------------------------------------------------------- *
      real function f_ntr(val_x,val_y,val_z)
      use m_usr
      implicit none
      !include 'usr.com'
      
      real val_x, val_y, val_z
      real norm
      real eps, delta, t
*     MdT  Adapt eps to system with dimensionless coordinates!      
!     eps   = par(SPL2)
      eps   = (r0dim/hdim)*par(SPL2)
      norm  = sqrt(val_x**2+val_y**2)/eps
      delta =0.05*norm
*     MdT  New function definition, to filter out steep slopes!
!     if ((val_z.ge.0.0).or.(val_z.lt.-norm)) then
!     f_ntr=0.0
!     elseif ((val_z.ge.-delta).and.(val_z.lt.0.0)) then
!     t=(val_z+delta)/delta
!     f_ntr=1.0-3.0*t**2+2.0*t**3
!     elseif ((val_z.ge.-norm).and.(val_z.lt.-norm+delta)) then
!     t=(val_z+norm)/delta
!     f_ntr=3.0*t**2-2.0*t**3
!     else
!     f_ntr=1.0
!     endif
      if (val_z.lt.-norm-delta) then
         f_ntr=1.0
      elseif ((val_z.ge.-norm-delta).and.(val_z.lt.-norm)) then
         t=(val_z+norm+delta)/delta
         f_ntr=1.0-3.0*t**2+2.0*t**3
      else
         f_ntr=0.0
      endif
      
      end
*     ---------------------------------------------------------------------------- *
      real function kd(x,x0)
      implicit none
      
      integer x,x0
      real val
      
      if (x.eq.x0) then
         kd=1.0
      else	
         kd=0.0
      endif
      end
*     ============================================================================ *
