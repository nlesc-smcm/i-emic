function [A] = inverseblockdiagonal(M, bs)
% M:  matrix
% bs: blocksize

  [n,m] = size(M);
  A = sparse([],[],[],n,m,0);

  tic
  fprintf('Creating block inverse... ');
  for i = 1:bs:n
	r = i:i+bs-1; %range
	A(r,r) = inv(M(r,r));
  end
  fprintf ('done (%f)\n', toc);

end
