%% For this to work you need to have libdap and loadapp and loadapp.m
%% in the path.

latmin = 10;
latmax = 80;
lonmin = 290;
lonmax = 355;
depth  = 4000;

n = 40; m = 40; l = 8;
qz = 1; % Not sure what to do with qz (yet)
periodic = false;

region    = 'natl'
mask_name = sprintf('%s_%d_%d_%d_%d',region,lonmin,lonmax,latmin,latmax);

if ~exist('temp')
  loaddap('http://iridl.ldeo.columbia.edu/SOURCES/.LEVITUS94/.ANNUAL/.temp/dods')
end

gridx = linspace(lonmin+.5,lonmax+.5,n);
gridy = linspace(latmin,latmax,m);
gridz = linspace(0,depth,l);

T = temp.temp;
%T(isnan(T)) = -999999;

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
	  JAC = [(T(d_j,d_i,d_k)-T(d_j,d_i-1,d_k))/(X(d_i)-X(d_i-1)),...
			 (T(d_j,d_i,d_k)-T(d_j-1,d_i,d_k))/(Y(d_j)-Y(d_j-1)),...
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

%transform_mask(mask_name, periodic);

figure(1);
%contourf(gridx,gridy,Tgrid(:,:,1)',10);
pcolor(gridx,gridy,Tgrid(:,:,1)')
xlim([lonmin,lonmax]);
ylim([latmin,latmax]);


figure(2)
%contourf(X,Y,T(:,:,1),10);
dimin = min(d_i_arr);
dimax = max(d_i_arr);
djmin = min(d_j_arr);
djmax = max(d_j_arr);

pcolor(X,Y,T(:,:,1));
xlim([lonmin,lonmax]);
ylim([latmin,latmax]);

figure(3)
pcolor(gridx,gridy,mask(:,:,1));
xlim([lonmin,lonmax]);
ylim([latmin,latmax]);


