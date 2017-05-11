[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');
surfm      = landm(2:n+1,2:m+1,l+1);    %Only interior surface points
landm_int  = landm(2:n+1,2:m+1,2:l+1);

srf = [];
greyness = .85;
srf(:,:,1) = (1-greyness*(surfm'));
srf(:,:,2) = (1-greyness*(surfm'));
srf(:,:,3) = (1-greyness*(surfm'));

atmos_nun = 1;
atmos_l = 1;
state  = readhdf5('atmos_output.h5', atmos_nun, n, m, atmos_l);

T0  = 15.0;   %//! reference temperature
RtD = 180/pi;

Ta  = reshape(T0 + state,n,m);
Tz  = mean(Ta,1); % zonal mean

%%
figure(10)

img = Ta';
contourf(RtD*x,RtD*(y),img,20,'Visible','off'); hold on;
image(RtD*x,RtD*(y),srf,'AlphaData',.2);
c = contour(RtD*x,RtD*(y),img,15,'Visible', 'on','linewidth',1);
colorbar
caxis([min(min(Ta)),max(max(Ta))])
hold off
drawnow
title('Atmospheric temperature')
xlabel('Longitude')
ylabel('Latitude')
exportfig('atmosTemp.eps')


%%
figure(11)

img = (Ta-repmat(Tz,n,1))';
contourf(RtD*x,RtD*(y),img,20,'Visible','off'); hold on;
image(RtD*x,RtD*(y),srf,'AlphaData',.2);
c = contour(RtD*x,RtD*(y),img,15,'Visible', 'on','linewidth',1);
colorbar
caxis([min(min(img)),max(max(img))])
hold off
drawnow
title('Atmospheric temperature anomaly')
xlabel('Longitude')
ylabel('Latitude')
exportfig('atmosTemp.eps')


