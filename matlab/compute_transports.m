function [] = compute_transports(solfile, datafile, use_mouse)

  if nargin < 3
	use_mouse = true;
	if nargin < 2
	  datafile = 'fort.44';
	  if nargin < 1
		solfile = 'fort.3';
	  end
	end
  end
  

  [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = ...
  readfort44(datafile);

  N = max([n,m]);
  maskp     = -sum(landm(2:end-1,2:end-1,2:end-1)-1,3)';
  int_maskp = maskp; %monitor
  tmp_maskp = maskp; %temporary

  coordinates = [];
  l = max(max(maskp));

  %% - READ SOLUTION - -------------------------------------------------
  [lab icp par xl xlp det sig sol solup soleig] = ...
  readfort3(la, solfile);

  %% - EXTRACT SOLUTION COMPONENTS - -----------------------------------
  [u,v,w,p,T,S] = extractsol(sol);

  %% - DEFINE CONSTANTS - ----------------------------------------------

  udim  = 0.1;                 %[m/s]   Velocity scale
  r0dim = 6.4e6;               %[m]     Radius of Earth
  T0    = 15;                  %[deg C] Reference temperature
  S0    = 35;                  %[psu]   Reference salinity
  RtD   = 180/pi;              %[-]     Radians to degrees  

  %% - Obtain coordinates for transect
  coords = getcoordinates(use_mouse, tmp_maskp);
  
  path = getpath(coords, int_maskp);

  % - Display path
  for i = 1:size(path,1);
	int_maskp(path(i,2),path(i,1)) = l+4*path(i,3);
  end
  imagesc(RtD*x,RtD*(y),int_maskp); set(gca,'ydir','normal');
  caxis([0,2*l])
  
  % - Compute transport for path
  transport = compute_transport(path, u*udim, v*udim, ...
								[x;xmax]*cos(yv(2:m+1))'*r0dim, ...
								[y;ymax]*r0dim, zw*hdim);

  disp(transport)
  
  hold on
  PSIB = bstream(u*udim,zw*hdim,[y;ymax]*r0dim);
  img = PSIB(2:end,:)';
  contour(RtD*x,RtD*(y),img,40,'linewidth',.5);
  hold off
  
end

function [transport] = compute_transport(path, u, v, x, y, z)

  N   = size(path,1);
  vel = zeros(N,1);
  dy  = y(2)-y(1); % grid increment is fixed in the y-direction

  for i = 1:size(path,1)	  

	if path(i,3) == 1 % u-orientation
	  % integrate u in depth
	  for j = 1:size(z)-1
		vel(i) = vel(i) + u(path(i,1),path(i,2),j) * (z(j+1)-z(j));
	  end
	  % append y grid increment
	  vel(i) = vel(i) * dy;
      
    elseif path(i,3) == 2 % v-orientation
      % integrate v in depth
      for j = 1:size(z)-1
    	vel(i) = vel(i) + v(path(i,1),path(i,2),j) * (z(j+1)-z(j));
      end
      % append x grid increment
	  dx = x(path(i,1)+1,path(i,2))-x(path(i,1),path(i,2));
      vel(i) = vel(i) * dx;
    end

  end

  % integrate over path and convert to Sverdrups
  transport = sum(vel)/1.e6; 
end

function [path] = getpath(coords, mask)
  x1 = coords(1,1);
  y1 = coords(1,2);
  x2 = coords(2,1);
  y2 = coords(2,2);

  % starting point
  x = x1;
  y = y1;

  % difference with target
  dx = (x - x2) / max(abs(x - x2),1);
  dy = (y - y2) / max(abs(y - y2),1);

  % initialize list that stores the path and orientation
  % orientation: or = 1 : u-direction
  %              or = 2 : v-direction

  or = 2;
  if dx == 0
	or = 1;
  end

  if 1%(mask(y,x) > 0)
	path = [x,y,or];
  else
	path = [];
  end

  corner = false;
  while (dx ~= 0) || (dy ~= 0)

	% traverse the transect, first in the x-direction
	if dx ~= 0
	  switch dx
		case 1
		  x = x - 1;
		case -1 
		  x = x + 1;
	  end
	  dx = (x - x2) / max(abs(x - x2),1);
	  or = 2;
	  if (dx == 0) && (dy ~= 0) 
		corner = true;
	  end

	% y-direction if x is aligned with x2
	elseif dy ~= 0
		   
	  switch dy
		case 1
		  y = y - 1;
		case -1 
		  y = y + 1;
	  end
	  dy = (y - y2) / max(abs(y - y2),1);
	  or = 1;
	  corner = false;
	end
	
	if ~corner && 1%(mask(y,x) > 0)
	  path = [path;x,y,or];
	end
  end
end

function [coords] = getcoordinates(use_mouse, mask)
  if use_mouse

	reset_mask = mask;
	fprintf('Usage: \n')
	fprintf('  Left mouse:    set coordinate\n')
	fprintf('  Right mouse:   retry\n')
	fprintf('  Middle mouse:  abort\n')

	coords = [];
	mx = max(max(mask));
	while size(coords,1) < 2
		  
	  imagesc(mask); set(gca,'ydir','normal');
	  
	  [x,y,button] = ginput(1);
	  x = round(x);
	  y = round(y);
	  if (button == 1)
		mask(y,x) = mx+8;
		coords = [coords;x,y];
	  elseif (button == 3)
		continue;
	  elseif (button == 2)
		break;
	  end

	  if size(coords,1) == 2
		d = coords(2,:) - coords(1,:);
		if d(1) ~= 0 && d(2) ~= 0
		  fprintf('We only allow straight meridional and zonal transects!\n')
		  fprintf('Resetting...\n');
		  mask = reset_mask;
		  coords = [];
		end
	  end
	end
	
	
	imagesc(mask); set(gca,'ydir','normal');
	pause(.2);
  else
	coords = [5,9;5,16]; % test coords
  end

  assert(norm(size(coords)-[2,2])==0);
end

