function [A] = load_atmosphere_matrix(basename)

  beg  = importdata([basename,'.beg']);
  ico  = importdata([basename,'.ico']);
  jco  = importdata([basename,'.jco']);

  n   = numel(beg)-1;
  nnz = beg(end)-1;

  ivals = zeros(nnz,1);
  jvals = jco;
  vals  = ico;

  row = 1;
  idx = 1;
  while row <= n
	for k = beg(row):beg(row+1)-1
	  ivals(idx) = row;
	  idx = idx + 1;
	end
	row = row + 1;
  end
  A = sparse(ivals, jvals, vals, n, n);
end
