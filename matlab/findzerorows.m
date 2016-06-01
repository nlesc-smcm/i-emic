function r = findzerorows(M)
  % M: matrix
  [m,n] = size(M);

  r = [];
  tic
  fprintf('Finding zero rows... ');
  for i = 1:m
	if (sum(abs(M(i,:))) == 0)
	   r = [r, i];
	end
  end
  fprintf ('done (%f)\n', toc);
	
end
