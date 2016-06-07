if ~exist('A') || force
  A = mmread('jacobian.ocean');
end

[N,M] = size(A);

[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');

n = 120;
m = 54;
l = 12;

assert(N == n*m*l*nun);

% (u,v,w,p,T,S)
iu  = 1:6:N;
iv  = 2:6:N;
iw  = 3:6:N;
ipt = 4:6:N;
ip  = [];
rp  = [];

SWS = load('singrows');
fprintf('building p indices...\n')
tic;
for i = 1:numel(SWS)
  if SWS(i) == 0
	 ip = [ip, ipt(i)];
	 rp = [rp, i];
  end
end
fprintf('building p indices... done %d\n', toc)
iuv = sort([iu,iv]);

Auv  = A(iuv,iuv);
Guv  = A(iuv,ip);
Duv  = A(ip, iuv);

if ~exist('iAuv') || force
  iAuv = inverseblockdiagonal(Auv,2);
end

E  = -Duv * Guv;
dE = diag(E);
for i = 1:numel(dE)
  if (abs(dE(i)) > 3000)
	fprintf('%6d  ', i);
	fprintf('%1.3e ', sum(full(Duv(i,:))));
	fprintf('\n');
	SWS(rp(i)) = 6;
  end
end
		 
C   = -Duv * iAuv * Guv;

C11 = C(rp<n*m,rp<n*m);
E11 = E(rp<n*m,rp<n*m);

condest(C11)

figure(1)
plot(diag(E11))
sc    = sum(abs(E11),2);
E11sc = diag(1./sc)*E11;

figure(2)
plot(diag(E11sc));
condest(E11sc)

figure(3)

c = diag(C);
z = zeros(numel(c),1);
for i = 1:numel(c);
  if (sum(abs(c(i))) < 1e-10)
	z(i) = 1;
  end  
end

fprintf('----------------------\n');
for i = 1:numel(ip)
  %fprintf('%3d      %3d %3d\n', i, z(i), SWS(rp(i)));
  if z(i) == 1 && SWS(rp(i)) == 0
	SWS(rp(i)) = 3;
  end
end

SWS = reshape(SWS,n,m,l);

figure(4)
imagesc(squeeze(SWS(:,:,end))')
set(gca,'ydir','normal')
