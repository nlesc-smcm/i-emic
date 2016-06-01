lvl = 1;

if ~exist('A') || force
  A = mmread('jacobian.ocean');
end

[N,M] = size(A);

n = 8; m = 8; l = 4; nun = 6;

assert(N == n*m*l*nun);

% (u,v,w,p,T,S)
iu  = 1:6:N;
iv  = 2:6:N;
iw  = 3:6:N;
ip  = 4:6:N;

iuv = sort([iu,iv]);

Auv  = A(iuv,iuv);
Guv  = A(iuv,ip);
Duv  = A(ip, iuv);

if ~exist('r') || force
  r = findzerorows(Duv);
end


r'-1
