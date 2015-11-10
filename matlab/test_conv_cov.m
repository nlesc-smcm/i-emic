%% This script tests the convergence of the covariance matrices'
%% eigenvalues

obs = [125,250,500,1000];

first = 1000;
eigs  = [];
for i = 1:numel(obs)
  [Q,M,V,D,sol] = build_covmat(first:first+obs(i)-1,[],10);
  eigs = [eigs ; diag(D)'];
end

save('eigs.mat','eigs')

for i = 1:numel(obs)
	EV(i,:) = eigs(i,:) / sum(eigs(i,:));
end
