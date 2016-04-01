%% For this to work you need to have libdap and loadapp and loadapp.m
%% in the path.

latmin = -80;
latmax = 80;
lonmin = 0;
lonmax = 360;
depth  = 4000;

n = 64; m = 32; l = 12;
qz = 1; % Not sure what to do with qz (yet)
periodic = true;

flat = true; % flat bathymetry
 
region    = 'global'
mask_name = sprintf('%s_%d_%d_%d_%d_%d_%d_%d',...
					region,lonmin,lonmax,latmin,latmax,...
					n,m,l);

if ~exist('temp')
  loaddap('http://iridl.ldeo.columbia.edu/SOURCES/.LEVITUS94/.ANNUAL/.temp/dods')
end

gridx = linspace(lonmin+.5,lonmax+.5,n);
gridy = linspace(latmin,latmax,m);
gridz = linspace(0,depth,l);

T = temp.temp;

if (flat)
  for i = 1:size(T,3)
	T(:,:,i) = T(:,:,1);
  end
end
   
if depth > Z(end-1)
  fprintf('this might give trouble\n');
end

Tgrid = zeros(n,m,l);
mask  = ones(m,n,l);

d_k_arr = [];
d_i_arr = [];
d_j_arr = [];

for k = 1:l
  Zdf = abs(Z-gridz(k));
  d_k = find(Zdf == min(Zdf),1,'first');
  for j = 1:m
	Ydf = abs(Y-gridy(j));
	d_j = find(Ydf == min(Ydf),1,'first');
	for i = 1:n
	  Xdf = abs(X-gridx(i));
	  d_i = find(Xdf == min(Xdf),1,'first');

	  % Create Jacobian
	  if (d_i-1 == 0 && periodic)
		d_i_m = size(T,2);
	  else
		d_i_m = d_i - 1;
	  end
	  if (d_j-1 == 0 && periodic)
		d_j_m = size(T,1);
	  else
		d_j_m = d_j - 1;
	  end
	  
	  JAC = [(T(d_j,d_i,d_k)-T(d_j,d_i_m,d_k))/(X(d_i)-X(d_i_m)),...
			 (T(d_j,d_i,d_k)-T(d_j_m,d_i,d_k))/(Y(d_j)-Y(d_j_m)),...
			 (T(d_j,d_i,d_k+1)-T(d_j,d_i,d_k))/(Z(d_k+1)-Z(d_k))];

	  dX  = [X(d_i)-gridx(i); ...
			 Y(d_j)-gridy(j); ...
			 Z(d_k)-gridz(k)];
	  
	  Tgrid(i,j,k) = T(d_j,d_i,d_k) + JAC*dX;

	  if (isnan(Tgrid(i,j,k)))
		mask(j,i,k) = 0;
	  end	  
	end
  end
end

save([mask_name,'.mat'], 'mask');

if exist('stop_here') && stop_here
   return
end

smooth_mask(mask_name,2,12);

transform_mask(mask_name, periodic);

figure(1);
%contourf(gridx,gridy,Tgrid(:,:,1)',10);
imagesc(gridx,gridy,flipud(Tgrid(:,:,1)'))
xlim([lonmin,lonmax]);
ylim([latmin,latmax]);
title('interpolated T');


figure(2)
%contourf(X,Y,T(:,:,1),10);
dimin = min(d_i_arr);
dimax = max(d_i_arr);
djmin = min(d_j_arr);
djmax = max(d_j_arr);

imagesc(X,Y,flipud(T(:,:,1)));
xlim([lonmin,lonmax]);
ylim([latmin,latmax]);
title('original T');

M = load([mask_name, '.mat']);
mask = M.mask;

figure(3)
imagesc(gridx,gridy,flipud(mask(:,:,1)));
title('mask level 1')

figure(4)
imagesc(gridx,gridy,flipud(mask(:,:,5)));
title('mask level 5')
