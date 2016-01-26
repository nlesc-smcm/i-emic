Jnum = [];

for i = 0:191
	T = importdata(['Fcol',num2str(i),'.ocean']);
	Jnum = [Jnum, T.data(:,3)];
end

thresh = 1e-12;

Jnum(Jnum < thresh) = 0;
Jnum = sparse(Jnum);

Jnumr = reordering(Jnum,6);
