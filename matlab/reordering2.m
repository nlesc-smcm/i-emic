function [res,P] = reordering2(mat)
dof  = 6;
N    = size(mat,1);    
idx  = 1:dof:N;
idx2 = sort([idx, idx+1]);
list = [idx2, idx+2, idx+3, idx2+4];

P   = speye(N);
P   = P(list,:);
res = P*mat*P';
end
