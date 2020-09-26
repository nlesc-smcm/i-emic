*******************************************************************************
      PROGRAM integration
* Version 3.1: AN 21/04/05
* Version 3.0: AN 28/02/03
* Version 2.0: HO 17/07/00
* Version 1.0: HD 15/05/99
* 
* New in Version 3.1: adapted to new block Gauss-Seidel preconditioner.
*
* New in Version 3.0: (mainly in subroutine step)
* - Smart choice of teta in scheme
* - Corrector: Newton-chord instead of Newton-Raphson
*   fix Jacobian over several time steps
* - Tuning of accuracy in solution of systems
* - Error estimate for local truncation error
* - Adaptive step size control
*
* A point from a previous time integration is read (or a starting point)
* and integrated in time for nstp steps. We keep the I/O structure
* exactly the same as in the continuation code, even reading/writing 
* eigenvalues/eigenvectors. 
* 
      USE m_start
      implicit none
      include 'bag.com'
      include 'res.com'
!      include 'start.com'
      real    g05caf
      integer iad
      integer irs, nstp, i, lab, ipert, outp, idum, ico
      real    time, ppert, upert, tdum

*  irs : Label of restart point
*  nstp: Number of time steps
*  ds  : Time-step size in dimensionless quantities
*  outp: Output control parameter
*  ipert : (0/1..npar) (N/Y) random perturbation of initial field 
*  ppert : parameter perturbation
*  upert : solution perturbation

      open(11,file= 'input')
      read(11,*) irs, nstp, ds, maxds, outp, 
     >             ipert, ppert, upert, idum, idum, tdum, ico
      close(11)

      open(98,status='replace')
      close(98)
      open(7,status='replace')
      close(7)

      CALL allocate_mat  ! allocate matrixarrays
      ALLOCATE(rb2old(ndim),r2o(ndim),dfdmu(ndim))
      ALLOCATE(rh0(ndim),rh1(ndim),rh2(ndim))

      call init
      p0 = 0
      ures = 0
      iad = 0

      IF (irs.EQ.0) THEN
         call stpnt(u)
         icp   = 0 
         time  = 0.0
      ELSE
         call read_data(4,irs,icp)
         time  = det
      END IF
      IF (ipert.ne.0) THEN
         DO i=1,ndim
            u(i) = u(i) + upert*g05caf()
         ENDDO
         par(ipert) = par(ipert) + ppert  
      ENDIF 
      call prpar(ds,nstp,irs)

      cstep = 0
      DO i = 1,nstp
         call step(time,ico)
         icp = icp + 1
         xl  = time 
         det = time 
         call monitor(0,0,7)
      ENDDO

      det = time
      call write_data(3,irs)

      CALL deallocate_mat  ! deallocate matrixarrays
      DEALLOCATE(rb2old, r2o, dfdmu, rh0,rh1,rh2)

      open(11,file= 'input')
      write(11,'(2i5,2e,i5,a,i2)') irs, nstp, ds, maxds, outp, 
     >          '   0 0 0 0 0 0 ',ico
      close(11)
      
      END
*************************************************************************
      SUBROUTINE step(time,ico)
      USE m_start
      USE m_gsolve
      USE m_bgspars
      implicit none
      include 'bag.com'
      include 'mix.com'
!      include 'start.com'
*     IMPORT/EXPORT
      real     time
      integer  ico
*     LOCAL
      integer  k,itmx,i
      real     rh(ndim),rb(ndim),res,eps,kapp,newds
      real     un(ndim),tdti,Bu(ndim)
      real     teta, error, D(ndim)
      integer  id
      logical  chst
*     FUNCTION
      real     linrm
*     
**HO  -B(u^n+1 - u^n)/dt + teta * rhs(u^n+1) + (1-teta) * rhs(u^n) = 0
**HO  teta = 0,.5,1 (Forward Euler, Trapezoidal, Backward Euler)
**HO  in the formulation below teta = 0 is not possible
**AN  teta = .5 + kapp*dt  with 0<kapp*dt<.5 (trunc.error u^n+1 is O(dt^3))
**AN  NB eps is the desired accuracy of the global solution
*     below A is the Jacobian matrix of rhs(u).
*
*     itmx : maximum number of iterations in Newton(-chord) method (<10)
*     eps  : desired accuracy after unit time
*            NB accuracy after time interval T is equal to eps*T

      chst = .true.
      itmx = 10
*     eps*T is acuracy of solution after time T
      eps  = 1.0e-04
      res  = 1.0e+00

      IF (cstep.EQ.0) THEN
         fix = .false.
         call rhs(u,rh2)
      END IF

      IF (fix) THEN
         rh0 = rh1
         rh1 = rh2
      ELSE
         rh0 = 0.0
         rh1 = rh2
      END IF

      kapp = 1.0
      kapp = min(kapp,.5/ds)
      teta = 1.0 ! .5 + kapp*ds
      tdti = 1.0 / (teta * ds)
*     calculate: - rhs(u^n)
      vmix_fix=0 ! ATvS-Mix
      call rhs(u,rh)
*     fill only matrix B
      call fillcolB
*     calculate total residu at n: tdti*Bu^n +rhs(u^n)
      call matBvec(u,Bu)
      rh = -((1-teta)/teta) * rh + tdti * Bu
*     compute secant guess
      un   = u + ds*up
*
      write(99,*) ' Newton process '
      nid = 0

      k=0
      DO WHILE ((k.LT.itmx).AND.(res.gt.eps))
         k = k+1
*     compute new preconditioner or not?
         fix = .false.
         IF (fix) THEN
            id = 2
         ELSE
            id = 0
         END IF
*     set tolerances for solver
         tolabs = 1.0e-06 
         tolred = 1.0e-06
*     compute first iterate for: -rhs(u^n+1)
         call rhs(un,rb)
         call matBvec(un,Bu)
*     compute total residue
         rb = rb - rh + tdti * Bu
*     fill matrix A - 2B/dt
         call matrix(un,tdti)
         call g_solve(rb,det,id)
*     correct pressure field
         call cor_pres(rb)
         un = un + rb
         res = linrm(rb,ndim)
         write(99,999) k,res
         IF ((k.EQ.itmx).AND.(res.GT.eps)) THEN
               STOP ' in step: Newton process not converged'
         ELSE
            fix = .true.
         END IF
         nid = nid + 1
      ENDDO

*     det = 1.0
*     compute approximation of derivative
      up = (un-u)/ds
      u = un
      time = time + ds
      cstep = cstep+1
*
 999  format(9x,'it: ',i2,' Res :',e12.4)
      END
*******************************************************************************
      SUBROUTINE prpar(ds,nst,irs)
      implicit none
      include 'bag.com'
*     IMPORT/EXPORT
      integer  nst,isw,irs,ide
      real     ds
*     LOCAL
      integer  i
*
*     open(7,position='append')
      write(7,996)
      write(7,999) (par(i),i=1,npar)
      write(7,996)
      write(7,998) irs,nst,ds
      write(7,996)
      write(7,995) 'time','max(psi)','min(psi)','max(psi)_at',
     >                    'min(psi)_at','ave(psi)_at'
*     close(7)
 999  format('*',5e15.5)
 998  format('*',4x,'irs = ',i4,4x,' nstp = ',i4,
     &       5x,'ds = ',e12.4)
 996  format('*')
 995  format(6x,6a12)
      END
***************************************************************
      SUBROUTINE pcorrect(un)
*     Set pressure in point intcond equal to zero.
      implicit none
      include 'usr.com'
*     IMPORT/EXPORT
      real    un(ndim)
      real    pref
*     LOCAL
      integer i,j,k,row
*     EXTERNAL
      integer find_row2

      pref = un(rowintcon)
      write(99,*) 'reference pressure: ', pref

      do k = 1, l
      do j = 1, m
      do i = 1, n
         if (landm(i,j,k) .eq. OCEAN) then
           row = find_row2(i,j,k,PP)
           un(row) = un(row) - pref
         endif
      enddo
      enddo
      enddo

      end
***************************************************************
      subroutine cor_pres(un)
      IMPLICIT none
      INCLUDE 'usr.com'

      real un(ndim), s1(ndim), s2(ndim), x1(ndim), x2(ndim), un2(ndim)
      real dp1, dp2
      integer k,j,i, find_row2, row
      integer px,py,pz,pv
      
      s1=0.0
      s2=0.0
      do k = 1,l
         do j = 1,m
            do i = 1,n
               if (landm(i,j,k)==ocean) then
                  row=find_row2(i,j,k,PP)
                  if (mod(i+j,2)==0) then
                     s1(row) = 1.0
                  else
                     s2(row) = 1.0
                  endif
               endif
            enddo
         enddo
      enddo
      
      s1 = s1/sqrt(sum(s1))
      s2 = s2/sqrt(sum(s2))

      dp1= dot_product(un,s1)
      dp2= dot_product(un,s2)
      call matAvec(s1,x1)
      call matAvec(s2,x2)
      un2 = un - dp1*s1 - dp2*s2

      un=un2

      end
      





