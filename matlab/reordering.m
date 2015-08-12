function [res,P] = reordering(mat,dof)
  if nargin < 2
    dof = 6;
  end
  N    = size(mat,1);    
  idx  = 1:dof:N;
  list = [];

  for i = 0:dof-1
    list = [list idx+i];
  end

  P = speye(N);
  P = P(list,:);
  res = P*mat*P';
end
