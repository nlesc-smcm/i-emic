%% For this to work you need to have libdap and loadapp and loadapp.m
%% in the path.

latmin = 10;
latmax = 75;
lonmin = 280;
lonmax = 360;
depth  = 4000;

if ~exist('sal')
  loaddap('http://iridl.ldeo.columbia.edu/SOURCES/.LEVITUS94/.ANNUAL/.sal/dods')
end

if ~exist('temp')
  loaddap('http://iridl.ldeo.columbia.edu/SOURCES/.LEVITUS94/.ANNUAL/.temp/dods')
end

if ~exist('temperature')
  loaddap('http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NODC/.WOA05/.Grid-1x1/.Annual/.mn/.temperature/dods')
end

figure(1)
S = sal.sal;
S(S<-100) = NaN;
contourf(X,Y,squeeze(S(:,:,1)),100); colorbar
xlim([lonmin,lonmax])
ylim([latmin,latmax])
title('levitus salinity')

figure(2)
T = temp.temp;
T(T<-100) = NaN;
contourf(X,Y,squeeze(T(:,:,1)),50); colorbar
xlim([lonmin,lonmax])
ylim([latmin,latmax])
title('levitus temperature')

figure(3)
S = sal.sal;

count = squeeze(sum(~isnan(S(:,lonmin:lonmax,:)),2));

S(isnan(S)) = 0;
IH = squeeze(sum(S(:,lonmin:lonmax,:),2))'./count';
contourf(Y,Z,(IH),50); colorbar
xlim([latmin,latmax])
ylim([0,4000])
title('levitus isohalines')
set(gca, 'ydir', 'reverse')

figure(4)
T = temp.temp;

count = squeeze(sum(~isnan(T(:,lonmin:lonmax,:)),2));
isnan(T(:,lonmin:lonmax,:));

T(isnan(T)) = 0;
TH = squeeze(sum(T(:,lonmin:lonmax,:),2))'./count';
contourf(Y,Z,(TH),50); colorbar
xlim([latmin,latmax])
ylim([0,depth])
title('levitus isothermals')
set(gca, 'ydir', 'reverse')
