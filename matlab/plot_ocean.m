%---------------------------------------------------------------------
% PLOTTHCM - Mother script for plotting THCM output
%
%  Father is M. den Toom, who conceived it 06-11-08     
%
%  Modified by Erik -> t.e.mulder@uu.nl
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
dx         = (xu(n+1)-xu(1))/n;
dy         = (yv(m+1)-yv(1))/m;
dz         = (zw(l+1)-zw(1))/l;

[qz,dfzt,dfzw] = gridstretch(zw);

%% - READ SOLUTION - -------------------------------------------------

[lab icp par xl xlp det sig sol solup soleig] = ...
readfort3(la,'fort.3');

%% - EXTRACT SOLUTION COMPONENTS - -----------------------------------
[u,v,w,p,T,S] = ...
extractsol(sol);

%% - INTEGRATIONS - --------------------------------------------------

% Barotropic streamfunction
PSIB = bstream(u*udim,zw*hdim,[y;ymax]*r0dim);

% Overturning streamfunction
PSIG = mstream(v*udim,[x;xmax]*cos(yv(2:m+1))'*r0dim,zw*hdim);
PSIG = [zeros(m+1,1) PSIG zeros(m+1,1)];

%% - CHECK SALINITY - ------------------------------------------------

check      = checksal(S,x,y,dfzt);
vol        = sum(sum(1-surfm).*cos(y'))*dx*dy;
fprintf(1,'Average salinity deficiency of %12.8f psu.\n', -check/vol) 

%% - PLOT THE RESULTS - ----------------------------------------------

figure(2)
contourf(RtD*xu,RtD*y,PSIB',15);
colorbar
title('Barotropic Streamfunction');
xlabel('Longitude')
ylabel('Latitude')
exportfig('bstream.eps')

figure(4)
contourf(PSIG',15);
colorbar
title('Overturning Streamfunction')
xlabel('Latitude')
ylabel('z (m)')
exportfig('mstream.eps')

figure(5)
Tp = T(:,:,l); 
contourf(RtD*xu(1:end-1),RtD*y,T0+Tp',15);
colorbar

figure(6)
Tp2 = squeeze(mean(T,1)); 
contourf(RtD*yv(1:end-1),z*hdim,Tp2'+T0,15);
colorbar
title('Isothermals')
xlabel('Latitude')
ylabel('z (m)')
exportfig('isopycnals.eps')
