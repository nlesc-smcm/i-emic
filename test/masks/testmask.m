%This is for a global 2deg mask with 12 vertical layers

if ~exist('A')
  A = mmread('jacobian.ocean');
end

[N,M] = size(A);

n = 120; m = 54; l = 12; nun = 6;

assert(N == n*m*l*nun);

% (u,v,w,p,T,S)
iu  = 1:6:N;
iv  = 2:6:N;
iuv = sort([iu,iv]);

Auv  = A(iuv,iuv);
%Auvr = reordering(Auv, 2);
%spy(Auvr)

spy(Auv)

condest(Auv)

Z = abs(diag(Auv))<1e3;
plot(Z)
