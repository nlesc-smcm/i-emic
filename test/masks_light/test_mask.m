if ~exist('A') || force
  A = mmread('jacobian.ocean');
end

[N,M] = size(A);

[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = ...
readfort44('fort.44');

assert(N == n*m*l*nun);

% (u,v,w,p,T,S)
iu  = 1:6:N;
iv  = 2:6:N;
iw  = 3:6:N;
ip = 4:6:N;

ZR  = load('zerorows');
SWS = load('singrows');

iuv = sort([iu,iv]);

Auv  = A(iuv,iuv);
Guv  = A(iuv,ip);
Duv  = A(ip, iuv);

if ~exist('iAuv') || force
  iAuv = inverseblockdiagonal(Auv,2);
end

C   = -Duv * iAuv * Guv;

c = diag(C);
z = zeros(numel(c),1);
for i = 1:numel(c);
  if (sum(abs(c(i))) < 1e-3)
	z(i) = 1;
  end  
end
dc11 = c(1:m*n);
DC11 = reshape(dc11, m, n);
SA11 = full(sum(A(iv(1:m*n),:),2));
SA11 = reshape(SA11, m, n);

figure(1)
imagesc(DC11')
set(gca,'ydir','normal')

figure(2)
imagesc(SA11')
set(gca,'ydir','normal')

figure(3)
plot(c(1:m*n))

fprintf('----------------------\n');
for i = 1:numel(ip)
  fprintf('%3d      %3d %3d\n', i, z(i), SWS(i));
  if z(i) == 1 && SWS(i) == 0
	SWS(i) = 3;
  end
end

SWS = reshape(SWS,n,m,l);

figure(4)
imagesc(squeeze(SWS(:,:,end))')
set(gca,'ydir','normal')
