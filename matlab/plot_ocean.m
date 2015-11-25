%---------------------------------------------------------------------
% PLOTTHCM - Mother script for plotting THCM output
%
%  Father is M. den Toom, who conceived it 06-11-08     
%
%  Modified by Erik, 2015 -> t.e.mulder@uu.nl
%---------------------------------------------------------------------

fprintf(1,'----------------------------------------------\n')

%% - DEFINE CONSTANTS - ----------------------------------------------

udim  = 0.1;                 %[m/s]   Velocity scale
r0dim = 6.4e6;               %[m]     Radius of Earth
T0    = 15;                  %[deg C] Reference temperature
S0    = 35;                  %[psu]   Reference salinity
RtD   = 180/pi;              %[-]     Radians to degrees

%% - READ MASK - -----------------------------------------------------

[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');
surfm      = landm(2:n+1,2:m+1,l+1);  %Only interior surface points
landm_int  = landm(2:n+1,2:m+1,2:l+1);
dx         = (xu(n+1)-xu(1))/n;
dy         = (yv(m+1)-yv(1))/m;
dz         = (zw(l+1)-zw(1))/l;

srf = [];
greyness = .85;
srf(:,:,1) = (1-greyness*(surfm'));
srf(:,:,2) = (1-greyness*(surfm'));
srf(:,:,3) = (1-greyness*(surfm'));

[qz,dfzt,dfzw] = gridstretch(zw);

%% - READ SOLUTION - -------------------------------------------------

[lab icp par xl xlp det sig sol solup soleig] = ...
            readfort3(la,'fort.3');

%% - EXTRACT SOLUTION COMPONENTS - -----------------------------------
[u,v,w,p,T,S] = extractsol(sol);

%% - INTEGRATIONS - --------------------------------------------------
% Barotropic streamfunction
PSIB = bstream(u*udim,zw*hdim,[y;ymax]*r0dim);

% Overturning streamfunction
PSIG = mstream(v*udim,[x;xmax]*cos(yv(2:m+1))'*r0dim,zw*hdim);
PSIG = [zeros(m+1,1) PSIG];

%% - CHECK SALINITY - ------------------------------------------------
check = checksal(S,x,y,dfzt);
vol   = sum(sum(1-surfm).*cos(y'))*dx*dy;
fprintf(1,'Average salinity deficiency of %12.8f psu.\n', -check/vol) 

%% Create Temperature
% build longitudinal average over non-land cells
Tl = zeros(m,l);
for k = 1:l
    for j = 1:m
        count = 0;
        for i=1:n            
            if landm_int(i,j,k) == 0
                count = count + 1;
                Tl(j,k) = Tl(j,k) + T(i,j,k);
            end
        end
        Tl(j,k) = Tl(j,k) / count;
    end
end

%% - PLOT THE RESULTS - ----------------------------------------------
figure(1)
img = PSIB(2:end,:)';
minval = min(min(img));
maxval = max(max(img));
contourf(RtD*x,RtD*(y),img,20,'Visible', 'off'); hold on;
imagesc(RtD*x,RtD*(y),img,'AlphaData',.5); hold on
image(RtD*x,RtD*(y),srf,'AlphaData',.9); hold on
contour(RtD*x,RtD*(y),img,20,'Visible', 'on','linewidth',2); hold off;
colorbar
caxis([minval,maxval])
title('Barotropic Streamfunction');
xlabel('Longitude')
ylabel('Latitude')
exportfig('bstream.eps')

%%
figure(2)
contourf(RtD*([y;ymax+dy/2]-dy/2),zw*hdim',PSIG',14);
colorbar
title('MOC (Sv)')
xlabel('latitude')
ylabel('depth (m)')
exportfig('mstream.eps',10,[20,7])

%%
figure(3)
Tp = T(:,:,l);
temp = flipud(T0 + Tp');
contourf(RtD*x,RtD*y,T0+Tp',15); hold on
% imagesc(RtD*x,RtD*y,temp); hold on
colorbar
contour(RtD*x,RtD*y,T0+1e-4*(surfm'),1,'k-','linewidth',2); hold off
title('Surface Temperature');
xlabel('Longitude');
ylabel('Latitude');

%%
figure(4)
Tp2 = squeeze(mean(T,1)); 
contourf(RtD*yv(1:end-1),z*hdim,Tl'+T0,15);
colorbar
title('Isothermals')
xlabel('Latitude')
ylabel('z (m)')
exportfig('isopycnals.eps')

%%
figure(5)
imagesc(flipud(PSIG')); 
grid on
title('Overturning Streamfunction')

%%
figure(6)
Sp = squeeze(mean(S,1)); 
contourf(RtD*yv(1:end-1),z*hdim,Sp'+S0,15);
colorbar
title('Isohalines')
xlabel('Latitude')
ylabel('z (m)')
exportfig('isopycnals.eps')

%%
figure(7)
Sp = S(:,:,l);
contourf(RtD*x,RtD*y,S0+Sp',15);
title('Surface salinity');
xlabel('Longitude');
ylabel('Latitude');
