function [out] = create_matrix_from_banded_storage(BA, dim, ksub, ksup)
  % Create dense matrix from banded storage BA
  % dim:  dimension of the matrix
  % ksub: number of lower diagonals
  % ksup: number of upper diagonals
		 
  ldba  = 2*ksub+ksup+1;     % leading dimension of banded storage
  out   = zeros(dim, dim);   % dense storage
  kdiag = ksub + 1 + ksup;   

  for icol = 1:dim
    i1 = max(1, icol-ksup);
    i2 = min(dim, icol+ksub);
    for irow = i1:i2
      irowb = irow - icol + kdiag;
      out(irow, icol) = BA(irowb, icol);
    end
  end
  
end
