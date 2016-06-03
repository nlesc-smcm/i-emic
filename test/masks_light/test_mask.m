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

if ~exist('iAuv') || force
  iAuv = inverseblockdiagonal(Auv,2);
end
C = -Duv * iAuv * Guv;

c = diag(C);
z = zeros(numel(c),1);
for i = 1:numel(c);
  if (sum(abs(c(i))) == 0)
	z(i) = 1;
  end  
end

SWS = load('singrows');
fprintf('----------------------\n');
for i = 1:numel(SWS)
  fprintf('%d      %d %d\n', i, z(i), SWS(i));
end
SWS = reshape(SWS,n,m,l);

figure(4)
imagesc(squeeze(SWS(:,:,end-1))')
set(gca,'ydir','normal')


