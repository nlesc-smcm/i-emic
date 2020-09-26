      integer ires
      real    ures, p0
      common /TRY/ ures(ndim), p0, ires
