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
ip  = 4:6:N;

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

condest(C)

figure(1)
plot(diag(C))

c = diag(C);
z = zeros(numel(c),1);
for i = 1:numel(c);
  if (sum(abs(c(i))) < 1e-3)
	z(i) = 1;
  end  
end

rg = m*n*l;


Prows = A(ip(1:rg),:);
Urows = A(iu(1:rg),:);
Vrows = A(iu(1:rg),:);
for i = 1:m*n
  idP = find(Prows(i,:));
  idU = find(Urows(i,:));
  idV = find(Urows(i,:));
end

Prows(44,:)
Prows(35,:)

dc11 = c(1:m*n);
DC11 = reshape(dc11, m, n);
SV11 = reshape(full(sum(A(iv(1:m*n),:),2)),n,m);
SU11 = reshape(full(sum(A(iu(1:m*n),:),2)),n,m);
SW11 = reshape(full(sum(A(iw(1:m*n),:),2)),n,m);
SP11 = reshape(full(sum(A(ip(1:m*n),:),2)),n,m);

figure(2)
imagesc(DC11')
title('Schur compl')
set(gca,'ydir','normal'); colorbar;

figure(3)
imagesc(SU11')
title('U')
set(gca,'ydir','normal'); colorbar;

figure(4)
imagesc(SV11')
title('V')
set(gca,'ydir','normal'); colorbar;

figure(5)
imagesc(SW11')
title('W')
set(gca,'ydir','normal'); colorbar;

figure(6)
imagesc(SP11')
title('P')
set(gca,'ydir','normal'); colorbar;


fprintf('----------------------\n');
for i = 1:numel(ip)
%  fprintf('%3d      %3d %3d\n', i, z(i), SWS(i));
  if z(i) == 1 && SWS(i) == 0
	SWS(i) = 3;
  end
end

SWS = reshape(SWS,n,m,l);
