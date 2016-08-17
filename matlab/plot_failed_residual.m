[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');
surfm = landm(2:n+1,2:m+1,l+1);  %Only interior surface points

%resa = load('residual.atmos');
%reso = load('residual.ocean');

%resa = load('solution.atmos');
%reso = load('solution.ocean');

%resa = load('failed_rhs.atmos');
%reso = load('failed_rhs.ocean');
%reso = load('xdot.ocean');
reso = load('stateDir.ocean');

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

norm(reso)
%norm(resa)
%norm([reso;resa]')s

L = l;

Uo = zeros(m,n); % zonal velocity
Wo = zeros(m,n); % vertical velocity	
To = zeros(m,n); % temperature
TL = zeros(m,n); %
So = zeros(m,n); %
Ta = zeros(m,n); %
Po = zeros(m,n); %

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

dim = m*n*l*nun;

for j = 1:m
  for i = 1:n
	Uo(m-j+1,i) = reso(find_row(nun,n,m,l,i,j,l,1));	
	Wo(m-j+1,i) = reso(find_row(nun,n,m,l,i,j,l,3));	
	To(m-j+1,i) = reso(find_row(nun,n,m,l,i,j,l,5));	
	TL(m-j+1,i) = reso(find_row(nun,n,m,l,i,j,L,5));
	LM(m-j+1,i) = surfm(i,j);
	So(m-j+1,i) = reso(find_row(nun,n,m,l,i,j,l,6));
	Po(m-j+1,i) = reso(find_row(nun,n,m,l,i,j,l-1,4));
  %	Ta(m-j+1,i) = resa(i + (j-1)*m);	
  end
end

mx = max(max(abs(To)));
[j,i] = find(abs(To) == mx);
fprintf('maximum in To at i=%d, j=%d\n', i, m-j+1);

figure(1);
v = [1 1];
imagesc(To); title('surface temperature residual'); colorbar
hold on;
contour(mx*P,[-1e-10,0,1e-10])
if norm(LM)
  contour(mx*LM,2);
end
hold off
exportfig('failed_residual_surftemp.eps',12,[20,12]);


figure(2);
mx = max(max(abs(TL)));
imagesc(TL); title(['temperature residual in layer ', num2str(L)]); colorbar
hold on;
contour(mx*P,[-1e-10,0,1e-10])
if norm(LM)
  contour(mx*LM,2);
end
hold off
exportfig('failed_residual_layer.eps',12,[20,12]);

figure(3); imagesc(So); title('surface salinity residual');  colorbar
mx = max(max(abs(So)));
hold on;
contour(mx*P,[-1e-10,0,1e-10])
if norm(LM)
  contour(mx*LM,2);
end
hold off
exportfig('failed_residual_salinity.eps',12,[20,12]);

%figure(4); imagesc(Ta); title('atmosphere temperature residual');colorbar
%mx = max(max(abs(Ta)));
%hold on;
%contour(mx*P,[-1e-10,0,1e-10])
if norm(LM)
%contour(mx*LM,2);
end
%hold off
%exportfig('failed_residual_atmos.eps',12,[20,12]);

figure(5); imagesc(Uo); title('surface zonal velocity');colorbar
mx = max(max(abs(Uo)));
hold on;
contour(mx*P,[-1e-10,0,1e-10])
if norm(LM)
  contour(mx*LM,2);
end
hold off
exportfig('failed_residual_zonal.eps',12,[20,12]);

figure(6); imagesc(Po); title('surface pressure');colorbar
mx = max(max(abs(Uo)));
hold on;
contour(mx*P,[-1e-10,0,1e-10])
if norm(LM)
  contour(mx*LM,2);
end
hold off

figure(7); imagesc(Wo); title('vertical velocity');colorbar
