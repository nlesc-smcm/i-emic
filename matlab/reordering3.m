function [res,P] = reordering3(mat)
	idx = symrcm(mat);
	N   = size(mat,1);
	P   = speye(N);
	P   = P(idx,:);
	res = P*mat*P';
end
