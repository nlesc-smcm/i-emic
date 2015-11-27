function [A,M] = load_ocean_matrix()
  %*************************************************
  % Function that loads the THCM output files into
  %  matrices A and M
  % Author: Erik Mulder -> t.e.mulder@uu.nl
  %*************************************************

  Abeg  = load('A.beg');
  Aco   = load('A.co');
  Ainfo = load('A.info');
  Ajco  = load('A.jco');
  %Arl   = importdata('A.rl');
  Bco   = load('B.co');

  n   = Ainfo(1);
  nnz = Ainfo(2);

  ivals = zeros(nnz,1);
  jvals = Ajco;
  vals  = Aco;

  row = 1;
  idx = 1;
  while row <= n
	for k = Abeg(row):Abeg(row+1)-1
	  ivals(idx) = row;
	  idx = idx + 1;
	end
	row = row + 1;
  end
  A = sparse(ivals, jvals, vals, n, n);
  M = sparse(1:n, 1:n, Bco, n, n);
end
