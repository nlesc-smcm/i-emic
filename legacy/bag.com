      include 'par.com'

      integer, parameter :: nf   =  4

      integer icp
      real    u, up, w, xl, xlp, det, sig, tval
      common /VEC/ u(ndim), up(ndim), w(ndim,nf)
     + , xl, xlp, det, sig(nf,2), tval, icp

