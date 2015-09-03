[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');
surfm      = landm(2:n+1,2:m+1,l+1);    %Only interior surface points
landm_int  = landm(2:n+1,2:m+1,2:l+1);

srf = [];
greyness = .85;
srf(:,:,1) = (1-greyness*(surfm'));
srf(:,:,2) = (1-greyness*(surfm'));
srf(:,:,3) = (1-greyness*(surfm'));

%otemp  = importdata('atmos_oceanTemp.txt');
state  = importdata('atmos_state.txt');
%sol   = importdata('atmos_sol.txt');

T0  = 15.0;   %//! reference temperature
RtD = 180/pi;

%To  = reshape(T0 + otemp,n,m);
Ta  = reshape(T0 + state,n,m);
%Ts  = reshape(T0 + sol,n,m);

% figure(5)
% contourf(RtD*x, RtD*y,To',15);
% colorbar
% title('SST ')
% xlabel('Longitude')
% ylabel('Latitude')
% exportfig('oceanTemp.eps')
%%
colormap jet
figure(6)
img = Ta';
contourf(RtD*x,RtD*(y),img,20,'Visible', 'off'); hold on;
imagesc(RtD*x,RtD*(y),img,'AlphaData',.8);
image(RtD*x,RtD*(y),srf,'AlphaData',.5); 
c = contour(RtD*x,RtD*(y),img,10,'Visible', 'on','linewidth',2); 
colorbar
caxis([min(min(Ta)),max(max(Ta))])
hold off


drawnow


title('Atmosphere')
xlabel('Longitude')
ylabel('Latitude')
exportfig('atmosTemp.eps')

% figure(7)
% contourf(RtD*x,RtD*y,Ts',15);
% colorbar
% title('Solution')
% xlabel('Longitude')
% ylabel('Latitude')
% exportfig('atmosTemp.eps')