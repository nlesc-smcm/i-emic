**************************************************
      SUBROUTINE grid
      implicit none
      include 'usr.com'
      integer i,j,k
*     EXTERNAL
      real  dfdz,fz

      write(99,*) '============GRID==========='
      write(99,10) xmin*180/pi,xmax*180/pi,ymin*180/pi,ymax*180/pi
 10   format(1x,'configuration:[',f6.1,1x,f6.1,'] x [',f6.1,1x,f6.1,']')
      write(99,20) n,m,l
 20   format(1x,'resolution:',i8,'x',i8,'x',i8)
      write(99,30) qz
 30   format(1x,'stretching parameter',g12.4)
      write(99,*) '============GRID==========='

      IF (abs(xmax-xmin-2*pi).LT.1.0e-3) periodic = .true.
      dx = (xmax-xmin)/N
      dy = (ymax-ymin)/M
      dz = (zmax-zmin)/L
      if (la.eq.1) then
         periodic =.true.
      endif
      DO i=1,n
        x(i) = (real(i)-0.5)*dx + xmin
        xu(i)= (real(i)    )*dx + xmin
      ENDDO
      xu(0) = xmin
      DO j=1,m
        y(j) = (real(j)-0.5)*dy + ymin
        yv(j)= (real(j)    )*dy + ymin
      ENDDO
      yv(0) = ymin
      DO k=1,l
        ze(k)  = (real(k)-0.5)*dz + zmin
        zwe(k) = (real(k)    )*dz + zmin
        z(k)   = fz(ze(k) ,qz)
        zw(k)  = fz(zwe(k),qz)
* compute derivatives of mapping at T- points
        dfzT(k) = dfdz(ze(k),qz)
* compute derivatives of mapping at w- points
        dfzW(k) = dfdz(zwe(k),qz)
      ENDDO
      zw(0) = zmin
      dfzw(0) = dfdz(zmin,qz)
      end
**************************************************
      REAL FUNCTION fz(z,qz)
      implicit none
* INPUT/OUTPUT
      real z,qz
* LOCAL
      real th,tth
      th = tanh (qz * (z+1))
      tth = tanh (qz)
      IF (qz.gt.1.0) THEN
      fz = -1 + th/tth
* for testing 
      ELSE
*      fz = z
      fz = z + (1.-qz) * z * (1-z)
      ENDIF
      RETURN
      END
************************************************
      REAL FUNCTION dfdz(z,qz)
      implicit none
* INPUT/OUTPUT
      real z,qz
* LOCAL
      real tth,ch
      ch = cosh (qz * (z+1))
      tth = tanh (qz)
      IF (qz.gt.1.0) THEN
      dfdz = qz  / (tth * ch * ch)
      ELSE
*      dfdz = 1.0
      dfdz = 1.0 + (1.-qz)*(1.-2.*z)
      ENDIF
      RETURN
      END

