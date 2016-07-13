%---------------------------------------------------------------------
% PLOTTHCM - Mother script for plotting THCM output
%
%  Father is M. den Toom, who conceived it 06-11-08     
%
%  Modified by Erik, 2015 -> t.e.mulder@uu.nl
%---------------------------------------------------------------------
atlantic = true;
glbl     = true;

fprintf(1,'----------------------------------------------\n')

%% - DEFINE CONSTANTS - ----------------------------------------------

udim  = 0.1;                 %[m/s]   Velocity scale
r0dim = 6.4e6;               %[m]     Radius of Earth
T0    = 15;                  %[deg C] Reference temperature
S0    = 35;                  %[psu]   Reference salinity
RtD   = 180/pi;              %[-]     Radians to degrees

%% - READ MASK - -----------------------------------------------------

[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = ...
readfort44('fort.44'); 

surfm      = landm(2:n+1,2:m+1,l+1);  %Only interior surface points
landm_int  = landm(2:n+1,2:m+1,2:l+1);
dx         = (xu(n+1)-xu(1))/n;
dy         = (yv(m+1)-yv(1))/m;
dz         = (zw(l+1)-zw(1))/l;

if n < 180  % hack
    atlantic = false;
	glbl     = false;
end    

if atlantic 
  % Atlantic horizontal range for 2deg landmask
  surfm_atl = imread('surfm_atl.bmp');
  atl_j = find(RtD*y > -40 & RtD*y < 60);
end

% change longitudinal view
div   = floor(n/4);
natl  = [div:n,1:div-1]; % NA focus
pacf  = 1:n;             % PA focus

if glbl
  range = natl;
else
  range = pacf;
end

srf = [];
greyness = .5;
srf(:,:,1) = (1-greyness*(interp2(surfm(range,:)','cubic')));
srf(:,:,2) = (1-greyness*(interp2(surfm(range,:)','cubic')));
srf(:,:,3) = (1-greyness*(interp2(surfm(range,:)','cubic')));
srf(srf<0) = 0;
srf(srf>1) = 1;

[qz,dfzt,dfzw] = gridstretch(zw);

%% - READ SOLUTION - -------------------------------------------------

[lab icp par xl xlp det sig sol solup soleig] = ...
readfort3(la,'fort.3');

%% - EXTRACT SOLUTION COMPONENTS - -----------------------------------
[u,v,w,p,T,S] = extractsol(sol);

if atlantic
  % Atlantic v
  v_atl = v;
  for j = 1:m
	for i = 1:n
	  if surfm_atl(i,j) ~= 0
		v_atl(i,j,:) = 0;
	  end
	end
  end
end

%% - INTEGRATIONS - --------------------------------------------------
% Barotropic streamfunction;
PSIB = bstream(u*udim,zw*hdim,[y;ymax]*r0dim);

% Overturning streamfunction
PSIG = mstream(v*udim,[x;xmax]*cos(yv(2:m+1))'*r0dim,zw*hdim);
PSIG = [zeros(m+1,1) PSIG];

if atlantic
  % Atlantic overturning streamfunction
  APSIG = mstream(v_atl(:,atl_j,:)*udim,[x;xmax]*cos(yv(atl_j))'*r0dim,...
				  zw*hdim);
  
  APSIG = [zeros(size(APSIG,1),1) APSIG];
end

%% - CHECK SALINITY - ------------------------------------------------
check = checksal(S,x,y,dfzt);
vol   = sum(sum(1-surfm).*cos(y'))*dx*dy;
fprintf(1,'Average salinity deficiency of %12.8f psu.\n', -check/vol) 

%% Create Temperature
% build longitudinal average over non-land cells
Tl = zeros(m,l);
Sl = zeros(m,l);
for k = 1:l
  for j = 1:m
    count = 0;
    for i=1:n            
      if landm_int(i,j,k) == 0
        count = count + 1;
        Tl(j,k) = Tl(j,k) + T(i,j,k);
		Sl(j,k) = Sl(j,k) + S(i,j,k);
      end
    end
    Tl(j,k) = Tl(j,k) / count;
	Sl(j,k) = Sl(j,k) / count;
  end
end

if atlantic
  %% Create Temperature
  % build longitudinal average over non-land cells
  Tla = zeros(m,l);
  for k = 1:l
	for j = 1:m
      count = 0;
      for i=1:n            
		if (landm_int(i,j,k) == 0) && (surfm_atl(i,j) == 0)
          count = count + 1;
          Tla(j,k) = Tla(j,k) + T(i,j,k);
		end
      end
      Tla(j,k) = Tla(j,k) / count;
	end
  end
end

%% - PLOT THE RESULTS - ----------------------------------------------
figure(1)
img = PSIB(2:end,:)';
minPSIB = min(PSIB(:));
maxPSIB = max(PSIB(:));
contourf(RtD*x,RtD*(y),img(:,range),20,'Visible', 'off'); hold on;
imagesc(RtD*x,RtD*(y),img(:,range),'AlphaData',.5); hold on
image(RtD*x,RtD*(y),srf,'AlphaData',.9); hold on
contours = linspace(minPSIB,maxPSIB,15);
contour(RtD*x,RtD*(y),img(:,range),contours,'Visible', 'on','linewidth',2); hold off;
colorbar
caxis([minPSIB,maxPSIB])
title('Barotropic Streamfunction');

if glbl
  xtl  = get(gca,'xticklabel');
  xtl2 = xtl;
  for i = 1:numel(xtl)
	xtl2{i} = num2str( str2num(xtl{i}) + round(RtD*x(div),-1) - 360 );
  end
  set(gca,'xticklabel',xtl2);
end

xlabel('Longitude')
ylabel('Latitude')
exportfig('bstream.eps',10,[50,25])

%%% 
figure(2)
contourf(RtD*([y;ymax+dy/2]-dy/2),zw*hdim',PSIG',40);
colorbar
title('MOC (Sv)')
xlabel('latitude')
ylabel('depth (m)')
exportfig('mstream.eps',10,[20,7])

%%
if atlantic
 figure(8)
 contourf(RtD*[y(atl_j);y(max(atl_j))+dy]-dy/2,zw*hdim',APSIG',14)
 colorbar
 title('AMOC (Sv)')
 xlabel('latitude')
 ylabel('depth (m)')
 exportfig('amstream.eps',10,[20,7])
end  

%% -------------------------------------------------------
figure(3)
Tsurf = T(:,:,l);
minT = T0+min(min(Tsurf));
maxT = T0+max(max(Tsurf));
%for j = 1:m
%  for i = 1:n
%	if surfm(i,j) == 1
%	  Tsurf(i,j) = -999;
%	end
%  end
%end

img  = T0 + Tsurf(range,:)';
%contourf(RtD*x,RtD*(y),img(:,range),20,'Visible', 'off'); hold on;
%set(gca,'color',[0.65,0.65,0.65]);
%image(RtD*x,RtD*(y),srf,'AlphaData',0.5); hold on
%contours = linspace(minT,maxT,20);
%contourf(RtD*x,RtD*(y),img,contours,'Visible', 'on','linewidth',1);
imagesc(RtD*x,RtD*(y),img)
set(gca,'ydir','normal');
%hold off

caxis([minT,maxT]);
colorbar
title('Surface Temperature', 'interpreter', 'none');
xlabel('Longitude');
ylabel('Latitude');

if glbl
  xtl  = get(gca,'xticklabel');
  xtl2 = xtl;
  for i = 1:numel(xtl)
	xtl2{i} = num2str( str2num(xtl{i}) + round(RtD*x(div),-1) - 360 );
  end
  set(gca,'xticklabel',xtl2);
end
exportfig('sst.eps',10,[50,25])

%% -------------------------------------------------------
figure(4)
contourf(RtD*yv(1:end-1),z*hdim,Tl'+T0,15);
%imagesc(RtD*yv(1:end-1),z*hdim,Tl'+T0);
set(gca,'ydir','normal');
%pcolor(RtD*yv(1:end-1),z*hdim,Tl'+T0);
colorbar
title('Temperature')
xlabel('Latitude')
ylabel('z (m)')
exportfig('isothermals.eps',10,[20,7])

if atlantic
  figure(5)
  contourf(RtD*yv(1:end-1),z*hdim,Tla'+T0,15);
  colorbar
  title('Atlantic temperature')
  xlabel('Latitude')
  ylabel('z (m)')
  exportfig('atl_isothermals.eps',10,[20,7])
end

%%
figure(6)
contourf(RtD*yv(1:end-1),z*hdim,Sl'+S0,15);
%imagesc(Sp'+S0);
set(gca,'ydir','normal')
%pcolor(RtD*yv(1:end-1),z*hdim,Sp'+S0);
colorbar
title('Isohalines')
xlabel('Latitude')
ylabel('z (m)')

exportfig('isohalines.eps',10,[20,7])

%%-----------------------------------------------------------------------------
figure(7)
Ssurf = S(:,:,l);
minS = S0+min(min(Ssurf));
maxS = S0+max(max(Ssurf));
%for j = 1:m
%  for i = 1:n
%	if surfm(i,j) == 1
%	   Ssurf(i,j) = -999;
%	end
%  end
%end

img  = S0 + Ssurf(range,:)';
% contourf(RtD*x,RtD*(y),img(:,range),20,'Visible', 'off'); hold on;
% set(gca,'color',[0.65,0.65,0.65]);
% image(RtD*x,RtD*(y),srf,'AlphaData',0.5); hold on
% contours = linspace(minS,maxS,40);
% contourf(RtD*x,RtD*(y),img,contours,'Visible', 'on','linewidth',1); hold off

imagesc(RtD*x,RtD*(y),img);
set(gca,'ydir','normal');
grid on
caxis([minS,maxS]);
colorbar
title('Surface Salinity', 'interpreter', 'none');
xlabel('Longitude');
ylabel('Latitude');

if glbl
  xtl  = get(gca,'xticklabel');
  xtl2 = xtl;
  for i = 1:numel(xtl)
	xtl2{i} = num2str( str2num(xtl{i}) + round(RtD*x(div),-1) - 360 );
  end
  set(gca,'xticklabel',xtl2);
end

exportfig('sss.eps',10,[50,25])
