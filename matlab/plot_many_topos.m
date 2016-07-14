% decide which plot to show in movie
str = ['select plottype: \n(0) surftemp, (1) barstreamfunc,',...
	   '(2) isothermals, (3) moc, (4) amoc (5) zonal surface vel.\n'];
plottype = input(str);

surftemp      = false;
barstreamfunc = false;
isothermals   = false;
moc           = false;
amoc          = false; % only working with 2deg landmask
zonalsurf     = false;

switch (plottype)
  case 0
	surftemp = true;
  case 1
	barstreamfunc = true;
  case 2
	isothermals = true;
  case 3
	moc = true;
  case 4
	amoc = true;
  case 5
	zonalsurf = true;
end

% Create array of strings with filenames of the states (UNIX)
[s,statenames] = system('ls -rt state_topo_*[0-9]* | sed "s/ / /" | sort ')
newlines = find(statenames == char(10));
filenames = [];
pars = [];
begin = 1;
k = 1;
for i = 1:numel(statenames)
  if i == newlines(k)
	filenames{k} = sprintf('%s',statenames(begin:i-1));
	pars(k) = str2num(sprintf('%s',statenames(begin+11:i-1)));
	k = k+1;
	begin = i + 1;
  end
end

filenames'
pars'

[s,masknames]  = system('ls -rt mask_*[0-9]* | sed "s/ / /" | sort ')
newlines = find(masknames == char(10));
maskfiles = [];
begin = 1;
k = 1;
for i = 1:numel(masknames)
  if i == newlines(k)
	maskfiles{k} = sprintf('%s',masknames(begin:i-1));
	k = k+1;
	begin = i + 1;
  end
end

maskfiles'

input('OK?\n')
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


% change longitudinal view
div   = floor(n/4);
div   = 1;
natl  = [div:n,1:div-1]; % NA focus
pacf  = 1:n;             % PA focus
range = natl;

srf = [];
greyness = .5;
srf(:,:,1) = (1-greyness*(interp2(surfm(range,:)','linear')));
srf(:,:,2) = (1-greyness*(interp2(surfm(range,:)','linear')));
srf(:,:,3) = (1-greyness*(interp2(surfm(range,:)','linear')));

[qz,dfzt,dfzw] = gridstretch(zw);

%--------------------------------------------------------------------
[~,~,~,~,~,~,~,sol2,~,~] = readfort3(la,filenames{end});

%% - EXTRACT SOLUTION COMPONENTS - ----------------------------------
[u,v,w,p,T,S] = extractsol(sol2);

minT = min(T(:)) + T0
maxT = max(T(:)) + T0

minU = min(u(:)) * udim;
maxU = max(u(:)) * udim;

%% - INTEGRATIONS - -------------------------------------------------

% Barotropic streamfunction
PSIB = bstream(u*udim,zw*hdim,[y;ymax]*r0dim);
minPSIB = min(PSIB(:))
maxPSIB = max(PSIB(:))

% MOC
PSIG = mstream(v*udim,[x;xmax]*cos(yv(2:m+1))'*r0dim,zw*hdim);
minPSIG = min(PSIG(:));
maxPSIG = max(PSIG(:));

if amoc
  % Atlantic horizontal range for 2deg landmask
  surfm_atl = imread('surfm_atl.bmp');
  atl_j = find(RtD*y > -40 & RtD*y < 70);

  % Atlantic v
  v_atl = v;
  for j = 1:m
	for i = 1:n
	  if surfm_atl(i,j) ~= 0
		v_atl(i,j,:) = 0;
	  end
	end
  end
  % Atlantic overturning streamfunction
  APSIG = mstream(v_atl(:,atl_j,:)*udim,[x;xmax]*cos(yv(atl_j))'*r0dim,...
				  zw*hdim);
  
  APSIG = [zeros(size(APSIG,1),1) APSIG];
  minAPSIG = min(APSIG(:));
  maxAPSIG = max(APSIG(:));
end

% get parameter
par2 = pars(end);
pare = par2;
%--------------------------------------------------------------------

[~,~,~,~,~,~,~,sol1,~,~] = readfort3(la,filenames{1});
par1 = pars(1);
parb = par1;

if surftemp
  fname = ['surftemp',sprintf('%5.4f',par2),'.avi']
elseif barstreamfunc
  fname = ['barstreamfunc',sprintf('%5.4f',par2),'.avi']
elseif isothermals
  fname = ['isothermal',sprintf('%5.4f',par2),'.avi']
elseif moc
  fname = ['moc',sprintf('%5.4f',par2),'.avi']
elseif amoc
  fname = ['amoc',sprintf('%5.4f',par2),'.avi']
elseif zonalsurf
  fname = ['zonalsurf',sprintf('%5.4f',par2),'.avi']
else
  fname = 'movie.avi'
end

writerObj = VideoWriter(fname, 'Motion JPEG AVI');
writerObj.FrameRate = 15;
writerObj.Quality = 90;
open(writerObj);
seconds = input('duration movie in seconds\n');
frames = writerObj.FrameRate * seconds;
par_incr = (pare - parb) / frames;
fhandle = figure('units','pixels','position',[0,0,1280,720]);
set(gca,'position',[0.05 0.1 .92 0.85],'units','normalized');
set(gca,'color','w','fontsize',15);

sol2 = sol1;
par2 = par1;
srf1 = [];
srf2 = [];
for file = 2:numel(filenames)
  sol1 = sol2;
  par1 = par2;

  % - READ SOLUTION - -------------------------------------------------
  [~,~,~,~,~,~,~, sol2,~,~] = readfort3(la,filenames{file});
  par2 = pars(file);

  % - GET LANDMASK ----------------------------------------------------
  [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,landm] = ...
  readfort44(maskfiles{floor(par2)+1});
  surfm = landm(2:n+1,2:m+1,l+1);  %Only interior surface points
  srf1(:,:,1) = (1-greyness*(interp2(surfm(range,:)','linear')));
  srf1(:,:,2) = (1-greyness*(interp2(surfm(range,:)','linear')));
  srf1(:,:,3) = (1-greyness*(interp2(surfm(range,:)','linear')));

  [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,landm] = ...
  readfort44(maskfiles{min(floor(par2)+2,numel(maskfiles))});
  surfm = landm(2:n+1,2:m+1,l+1);  %Only interior surface points
  srf2(:,:,1) = (1-greyness*(interp2(surfm(range,:)','linear')));
  srf2(:,:,2) = (1-greyness*(interp2(surfm(range,:)','linear')));
  srf2(:,:,3) = (1-greyness*(interp2(surfm(range,:)','linear')));

  delta = par2 - floor(par2);
  srf = (1-delta)*srf1 + delta*srf2;
  
  for pr = parb:par_incr:pare

	if pr < par1 || pr > par2
	  continue;
	end
	% HOMOTOPY --------------------------------------------------------
	a = (pr-par1) / (par2-par1);
	sol = (1-a)*sol1 + a*sol2;

	
	%% - EXTRACT SOLUTION COMPONENTS - -----------------------------------
	[u,v,w,p,T,S] = extractsol(sol);

	%% - INTEGRATIONS - --------------------------------------------------
	% Barotropic streamfunction
	if barstreamfunc
	  PSIB = bstream(u*udim,zw*hdim,[y;ymax]*r0dim);
	end

	if moc
	  % Overturning streamfunction
	  PSIG = mstream(v*udim,[x;xmax]*cos(yv(2:m+1))'*r0dim,zw*hdim);
	  PSIG = [zeros(size(PSIG,1),1) PSIG];
	end

	if amoc
	  % Atlantic v
	  v_atl = v;
	  for j = 1:m
		for i = 1:n
		  if surfm_atl(i,j) ~= 0
			v_atl(i,j,:) = 0;
		  end
		end
	  end
	  % Atlantic overturning streamfunction
	  APSIG = mstream(v_atl(:,atl_j,:)*udim,[x;xmax]*cos(yv(atl_j))'*r0dim,...
					  zw*hdim);
	  
	  APSIG = [zeros(size(APSIG,1),1) APSIG];
	end

	if zonalsurf
	   USURF = u(:,:,end) * udim;
	end

	%% Create Temperature
	% build longitudinal average over non-land cells
	if isothermals
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
	end
	
	if surftemp %----------------------------------

	  Tsurf = T(:,:,l);
	  for j = 1:m
		for i = 1:n
		  if surfm(i,j) == 1
			Tsurf(i,j) = -999;
		  end
		end
	  end
	  img  = T0 + Tsurf(range,:)';
	  contourf(RtD*x,RtD*(y),img(:,range),20,'Visible', 'off'); hold on;
	  set(gca,'color',[0.65,0.65,0.65]);
	  image(RtD*x,RtD*(y),srf,'AlphaData',0.5); hold on
	  contours = linspace(minT,maxT,40);
	  contourf(RtD*x,RtD*(y),img,contours); hold off
	  
	  colorbar
	  caxis([minT,maxT]);
	  title(['Surface Temperature ', sprintf('%5.4f',pr)], 'interpreter', 'none');
	  xlabel('Longitude');
	  ylabel('Latitude');
	  
	  xtl  = get(gca,'xticklabel');
	  xtl2 = xtl;
	  for i = 1:numel(xtl)
		xtl2{i} = num2str( str2num(xtl{i}) + round(RtD*x(div),-1) - 360 );
	  end
	  set(gca,'xticklabel',xtl2);
	  
	elseif barstreamfunc %----------------------------------

	  img = PSIB(2:end,:)';
	  minval = minPSIB;
	  maxval = maxPSIB;

	  contourf(RtD*x,RtD*(y),img(:,range),20,'Visible', 'off'); hold on;
	  %imagesc(RtD*x,RtD*(y),img,'AlphaData',.5); hold on
	  image(RtD*x,RtD*(y),srf,'AlphaData',.9); hold on
	  %contours = [fliplr(1-logspace(0.1,log10(abs(minval)),30)),-1+logspace(-.1,log10(abs(maxval)),10)]
	  contours = linspace(minval,maxval,40);
	  %contour(RtD*x,RtD*y,(surfm(range,:)'),1,'k-','linewidth',1); hold on
	  contour(RtD*x,RtD*(y),img(:,range),contours,'Visible', ...
			  'on','linewidth',2); hold on

	  hold off
	  colorbar
	  caxis([minval,maxval])
	  xlabel('Longitude')
	  title(['Barotropic Streamfunction ', sprintf('%5.4f',pr)], 'interpreter', 'none');
	  ylabel('Latitude')

      if div > 1
		xtl  = get(gca,'xticklabel');
		xtl2 = xtl;
	  
		for i = 1:numel(xtl)
		  xtl2{i} = num2str( str2num(xtl{i}) + round(RtD*x(div),-1) - 360 );
		end
		set(gca,'xticklabel',xtl2);
	  end

	  elseif moc
	  contours = linspace(minPSIG,maxPSIG,40);
	  contourf(RtD*([y;ymax+dy/2]-dy/2),zw*hdim',PSIG',contours,'linewidth',2);
	  colorbar
	  title(['MOC (Sv) ', sprintf('%5.4f',pr)],'interpreter', 'none');
	  xlabel('latitude')
	  ylabel('depth (m)')
	  
	elseif amoc % ---------------------------------------------------
	  contours = linspace(minAPSIG,maxAPSIG,40);
	  contourf(RtD*[y(atl_j);y(max(atl_j))+dy]-dy/2,zw*hdim',APSIG',contours,'linewidth',2);
	  colorbar
	  caxis([minAPSIG,maxAPSIG]);
	  title(['AMOC (Sv) ', sprintf('%5.4f',pr)], 'interpreter', 'none');
	  xlabel('latitude')
	  ylabel('depth (m)')

	elseif zonalsurf % ---------------------------------------------------
	  contours = linspace(minU,maxU,40);
	  imagesc(RtD*x,RtD*y, USURF');
	  set(gca,'ydir','normal');
	  colorbar
	  caxis([minU,maxU]);
	  title(['Surface zonal velocity (m/s) ', sprintf('%5.4f',pr)], 'interpreter', 'none');
	  xlabel('lon')
	  ylabel('lat')		   
	end
	
	fprintf('parameter: %f %f %f %f\n', pr, par1, par_incr, par2);	
	%----------------------------------
	frame = getframe(gcf);
	writeVideo(writerObj, frame);
  end
end
close(writerObj);
