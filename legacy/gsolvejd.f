      subroutine g_solve(RL,det,id)

*     This routines solves the equation Ax=rl for non-symmetric
*     matrix A using MRILU

      USE m_mat

      implicit none
      include 'usr.com'
!      include 'mat.com'
      include 'start.com'
      real,   dimension(ndim*(nun*np+1)):: coC 
      integer,dimension(ndim*(nun*np+1)):: jcoC 
      integer,dimension(ndim+1)         :: begC
      save coC, jcoC, begC
      real RL(ndim), det, nres, timeg1, timeg2
      integer id

      include 'mrilu.com'
      call cpu_time(timeg1)
*     SET TOLERANCES
      tolabs = 1.0e-06
      tolred = 1.0e-07
      AbsTol = tolabs
      RedTol = tolred

*      write(99,*) 'entering solve'
*     Print Out unscaled matrices
*     call writematrhs(RL)
*     stop
      if (id.eq.0) then
*     FIRST SYSTEM->PRECOND
*     Start from scracht, Zero the initial solution
         if( nid == 0 ) then
*            write(99,*) 'Writing out matrix rhs'
*            call writematrhs(RL)
*            call scaling(ndim, nun, begA, jcoA, coA, rs, cs)
            call scaling_land(begA, jcoA, coA, rs, cs)
*            cs = (/1.,1.,100.,1.,100.,10./)
*            rs = (/dx,dy,dz,dz*.01,.01,.1/)
            rsa = (/1.,1.,1.,1.,.3,1./)
            csa = (/1.,1.,1.,1.,.3,1./)
         endif
         call rcscale_a(rl, rs, cs, rsa, csa)
         coC  = coA
         jcoC = jcoA
         begC = begA
         call dsprav (ndim, 0.0D0, XS, 1)
*     Initialise work space buffer in 'Ibuf' and 'Dbuf'
         call iniwsb (bufsiz, Ibuf, Dbuf, ier)
         if (ier .LT. 0) stop 'in gsolve:  Failure in  iniwsb'
*     Construct the preconditioner and store location in 'ixSAVE'
         call cmpprc (ndim, nun, begC, jcoC, coC, ixSAVE,
     &        Ibuf, Dbuf, ier)
         if (ier .LT. 0) stop 'in gsolve:  Failure in  cmpprc'
      else
         call scalesol_a(RL,rs,rsa)
         IF (mod(id,2).EQ.0) THEN
*          FIRST SYSTEM->NOPRECOND
           call dsprav(ndim,0.0D0,XS,1)
         ELSE
*          SECOND SYSTEM->NO PRECOND
           XS = rb2old
           call scalesol(XS,1/cs)
         END IF
      endif
*     Solve the linear system  A x = RL using the preconditioner 'ixSAVE'
      call solprc (ndim, nun, begC, jcoC, coC, ixSAVE,
     &     XS, RL,
     &     Ibuf, Dbuf, ier)
      if (ier .LT. 0) stop 'in gsolve:  Failure in  solprc'
      call scalesol_a(XS,cs,csa)
*     call writevector(ndim,XS)
*     Set flag to use OLD PRECONDITIONING and OLD MATRIX
      id = id + 1
      det = 0.0D0
      RL = XS
      call cpu_time(timeg2)
      write(99,*) 'TOTAL SOLVE TIME: ',timeg2 - timeg1
      write(99,*) 'after precon, id = ', id 
      end









