F  = importdata('F.ocean');
F  = F.data(:,3);
J  = mmread('J.ocean');
Jr = reordering(J,6);


Jnum = [];

for i = 0:191
	T = importdata(['Fcol',num2str(i),'.ocean']);
	Jnum = [Jnum, T.data(:,3)];
end

thresh = 2e-17;

Jnum(Jnum < thresh) = 0;
Jnum = sparse(Jnum);

Jnumr = reordering(Jnum,6);

fprintf('nnz J    %d\n', nnz(J));
fprintf('nnz Jnum %d\n', nnz(Jnum));
