[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');

res  = load(fname);

% pargrid line color
pcol = [.3 .3 .3];

% parallel grid dimensions (horizontal)
npN = 1;
npM = 1;
Ndim = floor(n / npN) * ones(npN,1);
Mdim = floor(m / npM) * ones(npM,1);

for i = 1:mod(n,npN)
  Ndim(i) = Ndim(i)+1;
end

for j = 1:mod(m,npM)
  Mdim(j) = Mdim(j)+1;
end

assert(sum(Ndim) == n)
  
% create parallel checkerboard
P  = zeros(m,n);

irange = 0;
for i=1:npN
  irange = irange(end)+1:irange(end)+Ndim(i);
  jrange = 0;  
  for j = 1:npM
	jrange = jrange(end)+1:jrange(end)+Mdim(j);
	P(jrange,irange) = mod(i+j,2)-1/2;
  end
end

sol = zeros(nun,n,m,l+la);
idx = 1;
for k = 1:l+la
  for j = 1:m
	for i = 1:n
      for XX = 1:nun
		sol(XX,i,j,k) = res(idx);
		idx = idx + 1;
      end
	end
  end
end

[u,v,w,p,T,S] = extractsol(sol);

level = l;
figure(1); imagesc(u(:,:,level)'); title('u'); set(gca,'ydir','normal'); colorbar;
figure(2); imagesc(v(:,:,level)'); title('v'); set(gca,'ydir','normal'); colorbar;
figure(3); imagesc(w(:,:,level)'); title('w'); set(gca,'ydir','normal'); colorbar;
figure(4); imagesc(p(:,:,level)'); title('p'); set(gca,'ydir','normal'); colorbar;
figure(5); imagesc(T(:,:,level)'); title('T'); set(gca,'ydir','normal'); colorbar;
figure(6); imagesc(S(:,:,level)'); title('S'); set(gca,'ydir','normal'); colorbar;

landm_int = landm(2:end-1, 2:end-1, 2:end-1);

%figure(7); imagesc(l-sum(landm_int, 3)'); title('mask'); set(gca,'ydir','normal'); colorbar;

%pl = p(:,:,level);

fprintf('|u| = %e\n',norm(u(:)));
fprintf('|v| = %e\n',norm(v(:)));
fprintf('|w| = %e\n',norm(w(:)));
fprintf('|p| = %e\n',norm(p(:)));
%fprintf('|pl|= %e\n',norm(pl(:)));
fprintf('|T| = %e\n',norm(T(:)));
fprintf('|S| = %e\n',norm(S(:)));
fprintf('---------\n');
fprintf('|X| = %e\n',norm(res(:)));

