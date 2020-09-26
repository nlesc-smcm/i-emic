      real vmix_time
      integer vmix_row, vmix_col, vmix_dim
      integer vmix_ngrp, vmix_mingrp, vmix_maxgrp
      integer vmix_ipntr, vmix_jpntr
      integer vmix_flag, vmix_temp, vmix_salt
      integer vmix_fix, vmix_out, vmix_diff
      common /vmix/ vmix_time,
     +              vmix_row(ndim*(nun*np+1)), ! could be taken smaller: ndim*(2*15+1)
     +              vmix_col(ndim*(nun*np+1)), ! could be taken smaller: ndim*(2*15+1)
     +              vmix_dim, vmix_ngrp(ndim),
     +              vmix_mingrp,vmix_maxgrp,
     +              vmix_ipntr(ndim+1), vmix_jpntr(ndim+1),
     +              vmix_flag, vmix_temp, vmix_salt,
     +              vmix_fix, vmix_out, vmix_diff
