[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');

RtD   = 180/pi;              %[-]     Radians to degrees

res  = load(fname);
%tol = 1e-7;
%res(res>tol)  =  1;
%res(res<-tol)  = -1;
% pargrid line color
pcol = [.3 .3 .3];

surfacemask = landm(2:end-1,2:end-1,l+1);

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

land = (logical(p == 0)-.5);
figure(4); imagesc(p(:,:,level)'); title('p');
set(gca,'ydir','normal');
colormap(parula)
hold on
mxp = max(max(max(p)));
%contour(mxp*land(:,:,level)',4)
colorbar
hold off


figure(1); imagesc(RtD*x,RtD*(y),u(:,:,level)');
title('u'); set(gca,'ydir','normal');
colormap(parula)
colorbar
xlabel('Longitude')
ylabel('Latitude'); 

figure(2); imagesc(RtD*x,RtD*(y),v(:,:,level)');
title('v'); set(gca,'ydir','normal');
colormap(parula)
colorbar
xlabel('Longitude')
ylabel('Latitude'); 

figure(3); imagesc(w(:,:,level)'); title('w'); set(gca,'ydir','normal');
colormap(parula)
colorbar

figure(5); imagesc(T(:,:,level)'); title('T'); set(gca,'ydir','normal'); hold on
mxT = max(max(max(T)));
%contour(mxT*land(:,:,level)',4);
hold off
colormap(parula)
colorbar
figure(6); imagesc(S(:,:,level)'); title('S'); set(gca,'ydir','normal'); hold on
mxS = max(max(max(S)));
%contour(mxS*land(:,:,level)',4);
hold off

colormap(parula)
colorbar

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
