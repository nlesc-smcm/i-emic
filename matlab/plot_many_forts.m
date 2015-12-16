
% Create array of strings with filenames of the states
[s,statenames] = system('ls state*[0-9]* | sed "s/ / /" ')
newlines = find(statenames == char(10));

filenames = [];
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

srf = [];
greyness = .85;
srf(:,:,1) = (1-greyness*(surfm'));
srf(:,:,2) = (1-greyness*(surfm'));
srf(:,:,3) = (1-greyness*(surfm'));

[qz,dfzt,dfzw] = gridstretch(zw);

writerObj = VideoWriter( ['movie.avi'],'Motion JPEG AVI');
writerObj.FrameRate = 15;
writerObj.Quality = 90;
open(writerObj);
for file = 1:numel(filenames)
  %% - READ SOLUTION - -------------------------------------------------

  [lab icp par xl xlp det sig sol solup soleig] = ...
  readfort3(la,filenames{file});

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

  Tp = T(:,:,l);
  temp = flipud(T0 + Tp');
  contourf(RtD*x,RtD*y,T0+Tp',40); hold on
  % imagesc(RtD*x,RtD*y,temp); hold on
  colorbar
  caxis([2,22]);
  contour(RtD*x,RtD*y,T0+1e-4*(surfm'),1,'k-','linewidth',2); hold off
  title(['Surface Temperature ', filenames{file}], 'interpreter', 'none');
  xlabel('Longitude');
  ylabel('Latitude');
  %----------------------------------
  frame = getframe(gcf);
  writeVideo(writerObj, frame);
end
close(writerObj);
