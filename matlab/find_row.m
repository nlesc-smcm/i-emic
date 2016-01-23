function [out] = find_row(nun,n,m,l,irange,jrange,krange,XXrange)
  N = numel(irange)*numel(jrange)*numel(krange)*numel(XXrange);
  out = zeros(N,1);
  idx = 1;
  for k = krange
	for j = jrange
	  for i = irange
		for XX = XXrange
		  out(idx) = nun*((k-1)*n*m+ n*(j-1) + i-1)+XX;
		  idx = idx + 1;
		end
	  end
	end
  end
end					 

