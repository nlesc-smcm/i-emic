%*************************************************
% Script that loads the THCM output files into
%  matrices A and B
% Author: Erik Mulder -> t.e.mulder@uu.nl
%*************************************************

Abeg  = importdata('A.beg');
Aco   = importdata('A.co');
Ainfo = importdata('A.info');
Ajco  = importdata('A.jco');
Arl   = importdata('A.rl');
Bco   = importdata('B.co');

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
B = sparse(1:n, 1:n, Bco, n, n);
