function [J, co, ico, beg] = load_numjac(basename)

  co  = load([basename,'_co']);
  ico = load([basename,'_ico']);
  beg = load([basename,'_beg']);
  
  M = max(ico);
  N = numel(beg)-1;

  J = zeros(M,N);
  for j = 1:N
	for i = beg(j)+1:beg(j+1)
	  J(ico(i)+1,j) = co(i);
	end
  end
  J = sparse(J);
end
