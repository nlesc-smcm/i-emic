% decide which plots to show in movie
surftemp = false;
barstreamfunc = true;

% Create array of strings with filenames of the states (UNIX)
[s,statenames] = system('ls ocean_state*[0-9]* | sed "s/ / /" ')
newlines = find(statenames == char(10));
filenames = [];
paridx = 10; % continuation parameter, this should be obtained from filename
begin = 1;
k = 1;
for i = 1:numel(statenames)
  if i == newlines(k)
	filenames{k} = sprintf('%s',statenames(begin:i-1));
	k = k+1;
	begin = i + 1;
  end
end

filenames
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
[~,~, par2,~,~,~,~,sol2,~,~] = readfort3(la,filenames{end});
%% - EXTRACT SOLUTION COMPONENTS - ----------------------------------
[u,v,w,p,T,S] = extractsol(sol2);
minT = min(T(:)) + T0
maxT = max(T(:)) + T0
%% - INTEGRATIONS - -------------------------------------------------
% Barotropic streamfunction
PSIB = bstream(u*udim,zw*hdim,[y;ymax]*r0dim);
minPSIB = min(PSIB(:))
maxPSIB = max(PSIB(:))
% get parameter
par2 = par2(paridx);
%--------------------------------------------------------------------

[~,~, par1,~,~,~,~,sol1,~,~] = readfort3(la,filenames{1});
par1 = par1(paridx);

if surftemp
  fname = ['surftemp',sprintf('%5.4f',par2),'.avi']
elseif barstreamfunc
  fname = ['barstreamfunc',sprintf('%5.4f',par2),'.avi']
else
  fname = 'movie.avi'
end

writerObj = VideoWriter(fname, 'Motion JPEG AVI');
writerObj.FrameRate = 15;
writerObj.Quality = 90;
open(writerObj);
frames = writerObj.FrameRate * 20;
par_incr = (par2 - par1) / frames;
fhandle = figure('units','pixels','position',[0,0,1920,1080]);
set(gca,'position',[0.05 0.1 .92 0.85],'units','normalized');
set(gca,'color','w','fontsize',15);


sol2 = sol1;
par2 = par1;
for file = 2:numel(filenames)
  sol1 = sol2;
  par1 = par2;
  
  % - READ SOLUTION - -------------------------------------------------
  tic
  [~,~, par2,~,~,~,~, sol2,~,~] = readfort3(la,filenames{file});
  par2 = par2(paridx);
  toc;

  for pr = par1:par_incr:par2
	% HOMOTOPY --------------------------------------------------------
	a = (pr-par1) / (par2-par1);
	sol = (1-a)*sol1 + a*sol2;

	tic
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
	fprintf('data building done\n');
	toc

	tic
	if surftemp %----------------------------------
	  Tp = T(range,:,l);
	  temp = flipud(T0 + Tp');
	  contours = linspace(minT,maxT,40);
	  contourf(RtD*x,RtD*y,T0+Tp',contours); hold on
	  % imagesc(RtD*x,RtD*y,temp); hold on
	  colorbar
	  caxis([minT,maxT]);
	  contour(RtD*x,RtD*y,T0+1e-4*(surfm(range,:)'),1,'k-','linewidth',2); hold off
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
	  title(['Barotropic Streamfunction ', sprintf('%5.4f',pr)], 'interpreter', 'none');
	  xlabel('Longitude')
	  ylabel('Latitude')

	  xtl  = get(gca,'xticklabel');
	  xtl2 = xtl;
	  for i = 1:numel(xtl)
		  xtl2{i} = num2str( str2num(xtl{i}) + round(RtD*x(div),-1) - 360 );
	  end
	  set(gca,'xticklabel',xtl2);

	end
	fprintf('plotting done\n');
	fprintf('parameter: %f\n', pr);
	toc
	
	tic
	%----------------------------------
	frame = getframe(gcf);
	writeVideo(writerObj, frame);
	fprintf('frame export done\n')
	toc
  end
end
close(writerObj);
