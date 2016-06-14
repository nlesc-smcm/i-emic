A = mmread('jacobian.ocean');

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

iAuv = inverseblockdiagonal(Auv,2);
cond(full(iAuv))

rgs = find(~SWS);

C   = -Duv(rgs,:) * iAuv * Guv(:,rgs);
O   = -Duv(rgs,:) * Guv(:,rgs);

fprintf('cond iAuv              = %2.2e\n', cond(full(iAuv)));
fprintf('cond -Duv * iAuv * Guv = %2.2e\n', cond(full(C)));
fprintf('cond -Duv * Guv        = %2.2e\n', cond(full(O)));
fprintf('  min(diag(C))    = %2.2e\n', min(abs(diag(full(C)))));
fprintf('  max(diag(C))    = %2.2e\n', max(abs(diag(full(C)))));
fprintf('  max(C)          = %2.2e\n\n', max(max(abs(full(C)))));

fprintf('  min(Duv)        = %2.2e\n', min(min(abs(full(Duv(abs(Duv)>0))))));
fprintf('  max(Duv)        = %2.2e\n\n', max(max(abs(full(Duv(abs(Duv)>0))))));

fprintf('  min(diag(iAuv)) = %2.2e\n', min(abs(diag(full(iAuv)))));
fprintf('  min(iAuv)       = %2.2e\n', min(min(abs(full(iAuv(abs(iAuv)>0))))));
fprintf('  max(iAuv)       = %2.2e\n\n', max(max(abs(full(iAuv(abs(iAuv)>0))))));

fprintf('  min(Guv)        = %2.2e\n', min(min(abs(full(Guv(abs(Guv)>0))))));
fprintf('  max(Guv)        = %2.2e\n\n', max(max(abs(full(Guv(abs(Guv)>0))))));


figure(1)
plot(diag(C))

%vsm(C)
%vsm(O)

C   = -Duv * iAuv * Guv;
O   = -Duv * Guv;



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

DUV = reshape(full(sum(abs(Duv(1:m*n,:)),2)),n,m);
GUU = reshape(full(sum(abs(Guv(1:2:2*m*n,:)),2)),n,m);
GUV = reshape(full(sum(abs(Guv(2:2:2*m*n,:)),2)),n,m);

figure(7)
imagesc(DUV')
title('DUV')
set(gca,'ydir','normal'); colorbar;

figure(8)
imagesc(GUU')
title('GUU')
set(gca,'ydir','normal'); colorbar;

figure(9)
imagesc(GUV')
title('GUV')
set(gca,'ydir','normal'); colorbar;

SC2 = reshape(diag(C(1:m*n,:)),n,m);
SC3 = reshape(diag(C(1:m*n,:)),n,m);

figure(10)
imagesc(SC2')
title('SC2')
set(gca,'ydir','normal'); colorbar;

figure(11)
imagesc(SC3')
title('SC3')
set(gca,'ydir','normal'); colorbar;
