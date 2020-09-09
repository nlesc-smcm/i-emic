*******************************************************************************
      PROGRAM continu
      USE m_start
      implicit none
      include 'bag.com'
      include 'mix.com'
      include 'res.com'
!      include 'start.com'
      integer icpo, lab, irs, nstp, ana, ide, isw, isp, ise, i, file,ico
      real    ts

*  irs : Label of restart point
*  nstp: Number of steps along the branch
*  ds  : step size
*  maxds : maximal step size
*  icp : Bifurcation parameter
*  ana : Starting solution available (1/0)
*  ide : Detect special points       (1/0)
*  isw : Branch switching            (1/0)
*  isp : Eigenvalues during run      (1/0)
*  ise : Eigenvalues at end of run   (1/0)
*  tval: Target value for parameter icp (needed in userfun)
*  ico : Stable solution (postproc. only) (1/0)
*  ires: Residue continuation        (1/0/-1)

      call cpu_time(timea)
      open(11,file= 'input')
      read(11,*) irs,nstp,ds,maxds,icp,ana,ide,isw,isp,ise,tval,ico,ires
      close(11)
      if (ico .eq. 1) nstp = 0
      cstep=0
      call init
      IF (irs.EQ.0) THEN
         call stpnt(u)
         call initw(w,ndim,nf)
         icpo = 0
      ELSE
         call read_data(4,irs,icpo)
      END IF
      CALL allocate_mat  ! allocate matrixarrays
      ALLOCATE(rb2old(ndim),r2o(ndim),dfdmu(ndim))
      ALLOCATE(rh0(ndim),rh1(ndim),rh2(ndim))

! ATvS-Mix ---------------------------------------------------------------------              
c      call frs_test
! --------------------------------------------------------------------- ATvS-Mix

**
*  residue continuation 
**
      p0 = 0.0
      ures = 0.0
      select case(ires)
      case(0)
         call rhs(u,ures)
	 DO i=1,ndim
         write(15,16) ures(i)
	 ENDDO
 16      format(g14.7)
         STOP
      case(1)
         p0 = 1.0
         read(15,*) ures
      end select

      lab = irs
      par(ARCL)= 1.0/ndim
      xl = par(icp)
      call prpar(nstp,isw,ide,irs)
      IF ( (ana.NE.1).AND.(irs.EQ.0) ) call initialsol
      IF (icpo.NE.icp) call newdir
*
      call write_ind
      DO i = 1,nstp
         cstep=i
         call bstep(isw,ide,isp,lab)
         call outreg(u)
         IF (i.EQ.nstp) then
	   call write_data(3,lab)
	 else
	   rewind(3)
	   call write_data(3,lab)
	   lab=lab-1
	 endif  
c	 call write_frs(u)
      ENDDO

      IF (ise.EQ.1) call eigen

      if (nstp.EQ.0) then
      call monitor(0,1,7)
      call write_data(3,lab)
      endif

      CALL deallocate_mat ! deallocate matrixarrays
      DEALLOCATE(rb2old, r2o, dfdmu, rh0,rh1,rh2)

      open(11,file= 'input')
      write(11,'(2i5,2e,6i4,e,2i4)') lab,nstp,ds,maxds,icp
     >                         ,ana,ide,isw,isp,ise,tval,ico,ires
      close(11)
      call cpu_time(timeb)
      write(99,*) ' Newton total solve time:', timeb-timea
      write(99,'(a26,f10.3)') 'MIX| mixing time:         ', vmix_time
      write(99,'(a26,f10.3)') 'MIX| percent of total:    ',
     +  100.0*vmix_time/(timeb-timea)
      
      END
********************************************************************
      SUBROUTINE newdir
  
      USE m_gsolve
      USE m_bgspars
      USE m_start

*     Calculates a new branch derivative with respect to par(icp)
      implicit none
      include 'bag.com'
!      include 'start.com'
*     IMPORT/EXPORT
      real     rbp
*     LOCAL
      integer  i,id,ier
      real     rb(ndim,2),xnorm
*     FUNCTIONS
      real     l2nrm,zeta,dudl,duds,linrm,timem1,timem2
*
      write(99,*) 'entering newdir'
      zeta  = par(ARCL)
      call cpu_time(timem1)
      call rhsbrs(u,xl,rb,rbp)
      call matrix(u,0.0)
      call cpu_time(timem2)
      write(99,*) 'TIME TO CONSTRUCT MATRIX/RHS:           ',
     &           timem2 - timem1
      id    = 0
!      tolabs = 1.0e-05
!      tolred = 1.0e-06
      tolabs = 1.0e-08
      tolred = 1.0e-09
      call g_solve(rb(:,2),det,id)
*      call dcopy(ndim,rb(:,2),1,rb2old,1)
      rb2old = rb(:,2)
      dudl  = l2nrm(rb(:,2),ndim) 
      xnorm = sqrt( zeta*dudl*dudl+1.0)
      up    =-rb(:,2)/xnorm
      duds  = l2nrm(up,ndim) 
      xlp   = 1.0/xnorm
      write(99,999) duds,xlp
      write(99,*) 'Newdir done'
*
 999  format('up: ',e12.4,' lp  : ',e12.4)
      END
********************************************************
      SUBROUTINE bstep(isw,ide,isp,lab)
      USE m_start
      implicit none
      include 'bag.com'
      include 'mix.com'
!      include 'start.com'
*     IMPORT/EXPORT
      integer   isw,itp,ide,isp,lab
*     LOCAL
      real      deto,xlpo,levo,usro
      real      lev,usr,eps
      real      f0,f1,s0,s1
      real      dsw,dsold
      logical   syes,dyes,bifp
*     FUNCTIONS
      real      userfun
* 
      itp  = 0
      dyes = .false.
      syes = .false.
      eps  = 1.0e-06
      s0   = 0.0
      deto = det
      xlpo = xlp
      levo = sig(1,1)
      usro = userfun()
*
      call step(itp,lab)
      IF (isp.EQ.1) call eigen
      lev = sig(1,1)
      usr = userfun()
      s1  = s0 + ds
*
*     Check for singularity
*
      IF ( deto*det.LT.0.0 ) THEN
        write(99,*) ' Bifurcation point '
        itp=1
        f0 = deto
        f1 = det
        dyes=.true.
      ENDIF
      IF ( xlpo*xlp.LT.0.0 ) THEN
        write(99,*) ' Limit point '
        itp=2
        f0 = xlpo
        f1 = xlp
        dyes=.true.
      ENDIF
      IF ( levo*lev.LT.0.0 ) THEN
        IF ( abs(sig(1,2)).GT.eps ) THEN
          write(99,*) ' Hopf bifurcation point '
          itp=3
        ELSE
          write(99,*) ' Bifurcation point (ev) '
          itp=3
        ENDIF
        f0 = levo
        f1 = lev
        dyes=.true.
      ENDIF
      IF ( usro*usr.LT.0.0 ) THEN
        write(99,*) ' User defined point '
        itp=4
        f0 = usro
        f1 = usr
        dyes=.true.
      END IF
      IF (ide.EQ.0) dyes = .false.
      IF (dyes) call detect(s0,s1,f0,f1,itp,dyes)
*
      IF (isw.NE.0) syes = .true.
      IF ( (itp.EQ.1).OR.(itp.EQ.3) ) bifp = .true.
*     Branchswitch if isw>0, ipt=1 or itp=3 (bifp) and detection was
*     succesfull: dyes=.true.
      IF ( (syes).AND.(dyes).AND.(bifp) )THEN
        dsw=ds*isw
        call switch
        dsold=ds
        ds=dsw
        call step(itp,lab)
        dsw=ds
        ds=dsold
        call eigen
        itp=0
        ds=dsw
        isw=0 
      END IF
      call monitor(itp,0,7)
      IF ( (itp.eq.4).OR.(itp.eq.3) )THEN
        call write_data(3,lab)
        call cpu_time(timeb)
        write(99,*) ' Newton total solve time:', timeb-timea
        write(99,'(a26,f10.3)') 'MIX| mixing time:         ', vmix_time
        write(99,'(a26,f10.3)') 'MIX| percent of total:    ',
     +  100.0*vmix_time/(timeb-timea)
        STOP
      ENDIF
*
      END
********************************************************************
      SUBROUTINE switch
      implicit none
      include 'bag.com'
*     COMMON
      real     phl
      common /DD10/ phl(ndim)
*     LOCAL
      integer  i
      real     gmax,eps,phi1(ndim),rf1,rf2
*     FUNCTION
      real     linrm
*
      eps = 1.0e-06
      gmax = linrm(phl,ndim)
      write(99,999) gmax
      call kern(phi1)
      IF (gmax.LT.eps) THEN
       up = phi1
       xlp= 0.0
      ELSE
       rf1= 0.0
       rf2= 0.0
       DO i=1,ndim
         rf1=rf1+up(i)*phi1(i)
         rf2=rf2+up(i)*phl(i)
       ENDDO
       xlp =-rf1/(xlp-rf2)
       up = (phi1-xlp*phl)
      END IF
 999  format(' gmax: ',e12.4)
      END
***************************************************************
      SUBROUTINE step(itp,lab)

      USE m_gsolve
      USE m_bgspars
      USE m_start

*     Predicts solution and corrects it with Newton method
*     or adaptive Shamanskii method
      implicit none
      include 'bag.com'
      include 'mix.com'
!      include 'usr.com'
!     can not include both, because both include 'par.com'
*     IMPORT/EXPORT
      integer       itp,lab
*     COMMON
      real          phl
      common /DD10/ phl(ndim)
*     LOCAL
      integer       i,iter,itmx,id,nsolves,NCit,NEit
!      integer       Nopt,rest,NewtonMethod
      integer       rest
      real, DIMENSION(ndim,2) :: rb
      real                    :: rbp,dl,res,eps,Klimit
      real, DIMENSION(ndim)   :: un
      real      :: xln,time1,time2,timem1, timem2, tstart,tstop
      real          del2,del1,del0, dsold,xnorm,dudl
      logical       fini,gtKn,Newt
*     FUNCTIONS
      real          testfunc,l2nrm,zeta,duds
*     EXPLANATION OF IMPORTANT VARIABLES
*     itmx : maximum iterations
*     NCit : number of iteraions with same Jacobian
*     NEit : number of succesive Newton iterations
*     Nopt : optimal number of iterations (3<Nopt<11)
*     rest : counts succesive failing steps, restarts
*     fini : .true. if continuation step is finished
*     gtKn : if .true. value of Kn will be computed
*     meth : number of method that is used as corrector
*     Newt : .true. if Newton is used
*     res  : residual
*     del* : memory of residuals
*
      call cpu_time(time1)
      rest = 0
      fini = .false.
      call set_newton_pars
*     CHOOSE CORRECTOR METHOD in file 'usr.com'!
*     1:   Newton method
*     2:   Adaptive Shamanskii; Newton step if Newton-chord bad
*     3:   Newton-chord method
      IF (NM.EQ.1) THEN
         Newt = .true.
       ELSEIF (NM.EQ.2) THEN
         Newt = .true.
       ELSEIF (NM.EQ.3) THEN
         Newt = .false.
       END IF
      itmx = 2*Nopt+1
      zeta = par(ARCL)
*     
!      rb2old = 0
      DO WHILE (.NOT.fini)
         xln  = xl+ ds*xlp
         un   = u + ds*up
         par(icp) = xln
*     
         del2 = 1.0e+01
         del1 = 1.0e+02
         del0 = 1.0e+03
         eps  = 1.0e-03 ! ATvS
         res  = 1.0e+00
         NCit = 0
         NEit = 0
         iter = 0
         fix  = .false.
         gtKn = .false.
*     
         write(99,*) ' Newton process with step size',ds
         nid  = 0
         DO WHILE ((iter.LT.itmx).AND.(res.GT.eps))
            iter = iter+1
            IF (.NOT.fix) THEN
               NCit = 0
               NEit = NEit + 1
            ELSE
               NEit = 0
            END IF
            NCit = NCit+1
            call cpu_time(tstart)
*
            vmix_fix=0 ! ATvS-Mix
            CALL cpu_time(timem1)
            call rhsbrs(un,xln,rb,rbp)
            IF (.NOT.fix) THEN
               CALL matrix(un,0.0)
               id      = 0
               nsolves = 2
            ELSE
               id      = 2
               nsolves = 1
               rb(:,2) = rb2old
            END IF
            CALL cpu_time(timem2)
            WRITE(99,*) 'TIME TO CONSTRUCT MATRIX/RHS:           ',
     &           timem2 - timem1
*     
            DO i=1,nsolves
               IF (i.EQ.1) THEN
!                  tolabs = 1.0E-9
!                  tolred = 1.0E-10
                  tolabs = 1.0E-7
                  tolred = 1.0E-8
!                  tolabs = 1.0E-5
!                  tolred = 1.0E-6
               ELSE
!                  tolabs = 1.0E-8/del2
!                  tolred = 1.0E-9/del2
                  tolabs = 1.0E-7/del2
                  tolred = 1.0E-8/del2
!                  tolabs = 1.0E-4/del2
!                  tolred = 1.0E-5/del2
               END IF
               call g_solve(rb(1:ndim,i),det,id)
               IF ((.NOT.fix).AND.(i.EQ.2)) rb2old = rb(:,2)
            END DO
*     
            call solvdl(rb,rbp,un,xln,dl)
            call maxres(rb,dl,res,ndim)
            del0 = del1
            del1 = del2
            del2 = res
            Knc  = del2/del1
*      
            call cpu_time(tstop)
            IF ((NCit .EQ. 1).AND.(NM.NE.3)) THEN
               Tn = tstop-tstart
!	       gtKn = .true.
               IF (gtKn) THEN
                  Kn   = min(del2/(del1*del1),1D1)
                  gtKn = .false.
               END IF
               Klimit=2*Knc
            ELSEIF (NM.NE.3) THEN
               Tnc    = tstop-tstart
               IF (Kn.LE.0) Kn = 1
!               Klimit = exp(log(Kn*del2)*Tnc/Tn)
               ! ACdN: to prevent unnescesary Newton steps use: 
               Klimit = exp(log(max(Kn*del2,eps/del2))*Tnc/Tn)
            END IF
*     Check convergence and determine next step
            IF (NM.EQ.3) THEN
               write(99,999) iter,res,icp,xln,Knc
            ELSEIF (.NOT.fix) THEN
               write(99,998) iter,res,icp,xln,Kn
            ELSE
               write(99,997) iter,res,icp,xln,Knc,Klimit
            ENDIF
!            IF ((Knc.LT.1).OR.((Knc.LT.10).AND.(NCit.EQ.1))) THEN
            IF (Knc.LT.2) THEN
!               WRITE(99,*) "update solution "
               call update(un,rb,xln,dl,ndim)
               par(icp) = xln
            ELSEIF ((NM.EQ.2).AND.(NCit.GT.1)) THEN
               IF (Knc.LT.10) THEN
!                  WRITE(99,*) "update solution "
                  call update(un,rb,xln,dl,ndim)
                  par(icp) = xln      
               ELSE
                  Klimit = 1
                  del2 = del1
                  del1 = del0
                  iter = iter-1
               ENDIF
            ELSE
               iter = itmx
            END IF
            IF ((res.GT.1E1).AND.(NM.NE.2)) THEN
               iter = itmx
            END IF

            IF (NM.EQ.3) THEN
               fix = .true.
            ELSEIF (NM.EQ.1) THEN
               fix = .false.
            ELSEIF (Knc.GE.Klimit) THEN
               fix = .false.
            ELSEIF (((cstep.EQ.1).AND.(iter.EQ.1)).OR.(Kn.LE.1E-10)) THEN
               fix = .false.
            ELSE
               fix = .true.
            ENDIF
            IF ((.NOT.fix).AND.(NCit.EQ.1)) gtKn=.true.
*     
            nid = nid + 1
         END DO
         IF (res.LE.eps) THEN
            fini = .true.  
         ELSE IF (rest.GT.5) THEN
            write(99,*) '     Newton not converged after six restarts, s
     &aving previous'
            par(icp) = xl
            call write_data(3,lab)
            STOP
         ELSE IF ((itp.EQ.0).AND.(abs(ds).LT.1.0e-5)) THEN
            write(99,*) '     step size too small, saving previous'
            par(icp) = xl
            call write_data(3,lab)
            STOP
         ELSE
            write(99,*) '     Newton not converged: restart with smaller
     & step size' 
            rest = rest+1
         END IF
         dsold = ds
         IF ((itp.EQ.0).OR.(.NOT.fini)) call changestep(iter)
      END DO
*     update
      fix = .false.
      phl = rb(:,2)
      IF (PM.eq.1) THEN   ! Euler predictor
         dudl  = l2nrm(rb(:,2),ndim) 
         xnorm = sqrt( zeta*dudl*dudl+1.0)
         up    =-rb(:,2)/xnorm
         duds  = l2nrm(up,ndim) 
         xlp   = 1.0/xnorm
      ELSEIF (PM.eq.2) THEN ! Secant predictor
         up = (un-u)/dsold
         xlp= (xln-xl)/dsold
      ENDIF

      u = un
      xl= xln
      call cpu_time(time2)
      write(99,*) 'TOTAL E-N TIME: ',time2-time1
      duds = l2nrm(up,ndim)
      write(99,996) duds,xlp
*
 999  format(9x,'It:',i2,' Residu :',e12.4,'  Par(',i2,') :',e12.4,
     &'  Knc :',e12.4)
 998  format(9x,'It:',i2,' Residu :',e12.4,'  Par(',i2,') :',e12.4,
     &'  Kn  :',e12.4)
 997  format(9x,'It:',i2,' Residu :',e12.4,'  Par(',i2,') :',e12.4,
     &'  Knc :',e12.4,'  Klim :',e12.4)
 996  format('up  : ',e12.4,' lp  : ',e12.4)

      END
***************************************************************
      SUBROUTINE set_newton_pars
      USE m_start
      IMPLICIT none
      include 'usr.com'
      ! IMPORT/EXPORT
      NM = NewtonMethod
      Nopt = OptNewtonIter
      PM = predictor
      END
***************************************************************
      SUBROUTINE changestep(Niter)
*     Change step size ds depending on convergence
      USE m_start
      implicit none
*     IMPORT/EXPORT
      integer  Niter
*     LOCAL
      real     ksi,No,Ni
*
      No  = Nopt
      Ni  = Niter
      ksi = No/Ni
      IF (ksi.LT.5.0e-01) THEN
         ksi = 5.0e-01
      ELSE IF (ksi.GT.2) THEN
         ksi = 2
      END IF
      IF (abs(ksi*ds).LE.maxds) THEN
         ds = ksi*ds
      ELSE
         ds=maxds*(abs(ds)/ds)
      END IF
*
      END                 
***************************************************************
      SUBROUTINE initialz2(rb2)
*     Store initial value for z2
      USE m_start
      implicit none
      include 'bag.com'
!      include 'start.com'
*     IMPORT/EXPORT
      real rb2(ndim)
*
!      call dcopy(ndim,rb2,1,rb2old,1)
      rb2old = rb2
*
      END                 
***************************************************************      
      SUBROUTINE rhsbrs(un,xln,rb,rbp)
*     Calculate rhs-side for extended sytem (incl. arc-length par)
      USE m_start
      implicit none
      include 'bag.com'
!      include 'start.com'
*     IMPORT/EXPORT
      real     un(ndim),xln,rb(ndim,2),rbp
*     LOCAL
      integer  i
      real     delta,uu2,l2nrm,zeta
*
      zeta = par(ARCL)
      delta= 1.0e-06
      call rhs(un,rb(:,1))
*
*     arc-length parameterization:
      uu2 = 0.0
      DO i=1,ndim
*       uu2 = up(i)*(un(i)-u(i)) +uu2
        uu2 = (un(i)-u(i))*(un(i)-u(i)) +uu2
      ENDDO
      rbp = - uu2*zeta - (xln-xl)*(xln-xl) + ds*ds
*     rbp = - uu2*zeta - xlp*(xln-xl) + ds
*
*     Derivative of rhs with respect to bifurcation parameter
      if (fix .EQ. .false.) then
        par(icp) = par(icp)+delta
        call rhs(un,rb(:,2))
        par(icp) = par(icp)-delta
        DO i=1,ndim
          rb(i,2)=-(rb(i,2)-rb(i,1))/delta
        ENDDO
      endif
*
      END
***************************************************************
      SUBROUTINE solvdl(rb,rbp,un,xln,dl)
*     Solve for xln:
*     dl*2(xln-xl) = Rbp - du*2(un-u)
*     du = J^inv(rb1) - dl*J^inv(rb2)
      USE m_start
      implicit none
      include 'bag.com'
!      include 'start.com'
*     IMPORT/EXPORT
      real     rb(ndim,2),rbp,un(ndim),xln,dl
*     LOCAL
      integer  i
      real     inp1,inp2,zeta
*
      zeta = par(ARCL)
      inp1 = 0.0
      inp2 = 0.0         
      
      DO i=1,ndim
        inp1 = inp1 + rb(i,1)*2.0*(un(i)-u(i))*zeta
        inp2 = inp2 + rb(i,2)*2.0*(un(i)-u(i))*zeta
      ENDDO                            
      dl = (rbp-inp1)/(2.0*(xln-xl)-inp2)

*     UPDATE Delta u      
      DO i=1,ndim
        rb(i,1) = rb(i,1)-dl*rb(i,2)
      ENDDO
*
      END
********************************************************************
      SUBROUTINE update(un,rb,xln,dl,ndim)
*     update un and xln
      implicit none
*     IMPORT/EXPORT
      integer  ndim
      real     un(ndim),rb(ndim,2),xln,dl
*     LOCAL
      integer  i
*
      DO i=1,ndim
        un(i) = un(i)+rb(i,1)
      ENDDO
      xln = xln+dl
*
      END
********************************************************************
      SUBROUTINE maxres(du,dl,res,ndim)
      implicit none
*     IMPORT/EXPORT
      integer   ndim
      real      du(ndim),dl,res
*     LOCAL
      integer   i
      real      linrm
*
      res = linrm(du,ndim)
      IF (abs(dl).GT.res) res = abs(dl)
*
      END
******************************************************************
      SUBROUTINE detect(s0,s1,f0,f1,itp,dyes)
*     detect singular or user defined point using Secant technique
      USE m_start
      implicit none
      include 'bag.com'
!      include 'start.com'
*     IMPORT/EXPORT
      integer  itp
      real     s0,s1,f0,f1,f2
      logical  dyes
*     LOCAL
      integer  itmx,k
      real     eps,dsold
*     FUNCTION
      real     userfun
*    
      eps = 1.0e-06 ! ATvS-Mix
      itmx= 10       ! ATvS-Mix
*
      write(99,*) 'Detecting special point'
      DO k=1,itmx
         IF ( abs(f1).GT.eps ) THEN
*         IF ( abs(f1-f0).GT.eps ) THEN
            dsold = ds
            ds =-f1*(s1-s0)/(f1-f0)
            write(99,999) k,s1,f1
!            f0 = f1
!            s0 = s1
            call step(itp,0)
            IF (itp.EQ.1) f2 = det       ! Bifurcation point
            IF (itp.EQ.2) f2 = xlp       ! Limit point
            IF (itp.EQ.3) THEN           ! (Hopf) bifurcation (ev)
               call eigen
               f2 = sig(1,1)
            END IF
            IF (itp.EQ.4) f2 = userfun() ! User defined point
            
            IF (f2*f1.gt.0) THEN
               f1 = f2
               s1 = s1+ds
            ELSE
               f0 = f1
               f1 = f2
               s0 = s1
               s1 = s1+ds
            ENDIF
            ds = dsold
         else
            write(99,998) xl,f1
            dyes = .true.
            return
         ENDIF
      ENDDO
      write(99,*) 'no convergence in secant proces'
      dyes = .false.
*
 999  format(' Secant proces: It:',i2,' s: ',e12.4,' F(s): ',e12.4)
 998  format(' Special point at xl: ',e12.4,' F(xl): ',e12.4)
      END
********************************************************************
      SUBROUTINE kern(phi1)

      USE m_gsolve

      implicit none
      include 'bag.com'
*     IMPORT/EXPORT
      real    phi1(ndim)
*     LOCAL
      integer k,itmx,id
      real    phi0(ndim),res,eps,xnpt
*     FUNCTIONS
      real    l2nrm
*
      eps = 1.0e-06
      res = 1.0
      phi1= 1.0
      itmx= 5
*
      write(99,*) 'Inverse iteration: '
      call matrix(u,0.0)
      id = 0
      DO k=1,itmx
       IF (res.GT.eps) THEN
         phi0 = phi1
         call g_solve(phi1,det,id)
         xnpt = l2nrm(phi1,ndim)
         phi1 = phi1/xnpt
         phi0 = abs(phi1)-abs(phi0)
         res = l2nrm(phi0,ndim)
         write(99,999) res,xnpt
       ENDIF
      ENDDO
 999  format(1x,'Residu :',e12.4,2x,' Norm :',e12.4)
      END
********************************************************************
      real FUNCTION userfun()
      implicit none
      include 'bag.com'
      userfun = par(icp) - tval
      END
********************************************************************************
      SUBROUTINE prpar(nst,isw,ide,irs)
      USE m_start
      implicit none
      include 'bag.com'
!      include 'start.com'
*     IMPORT/EXPORT
      integer  nst,isw,irs,ide
*     LOCAL
      integer  i
*

      write(7,*) '---------------------------'
      write(7,'(5e15.5)') (par(i),i=1,npar)
      write(7,*) ' '
      write(7,998) irs,nst,isw,ide,ds
      write(7,*) ' '
!     write(7,'(40a)') 'typ par    val      Overturn(Sv)' 
!     write(7,*) 'itp par parval MOCG- MOCG+ MOCA- MOCA+ MOCP-
!    & MOCP+'
 998  format('    irs =',i3,', nstp =',i3,', isw =',i2,', ide =',i2,
     &       ',  ds = ',e12.4)

      end

********************************************************************
      SUBROUTINE initialsol
*     Calculates an initial solution if an analytical is not available    
      USE m_bgspars
      USE m_start
      implicit none
      include 'bag.com'
!      include 'start.com'
*     LOCAL
      real     dsold,pseudo

      up    = 0.0
      xlp   = 1.0
      dsold = ds
      ds    = 1.0e-08
      pseudo = par(ARCL)
      par(ARCL) = 0.0
      write(99,*) 'calculating initial solution'
      tolabs = 1.0E-6
      tolred = 1.0E-7
      call step(9,0)
      ds = dsold
      write(99,*) 'initial solution found'
      par(ARCL) = pseudo
      call monitor(9,0,7)

      end

***************************************************************
      SUBROUTINE updatez2(rb)
*     update z2 = dz2 + zold2
      USE m_start
      implicit none
      include 'bag.com'
!      include 'start.com'
*     IMPORT/EXPORT
      real rb(ndim)
*
!      call dcopy(ndim,rb,1,rb2old,1)
      rb2old = rb
*
      END
