function [transport, trpath] = compute_transports(solfile, datafile, trpath)

  transmode = 2;
  if nargin < 3
	transmode = 1;
	fprintf('Using mouse input to define transect/path\n');
	trpath = [];
  end
  if nargin < 2
	datafile = 'fort.44';
  end
  if nargin < 1
	solfile = 'fort.3';
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
  if transmode == 1
	coords = getcoordinates(tmp_maskp);  
	trpath = getpath(coords, int_maskp);
	
	%% - Display path with streamfunction
	for i = 1:size(trpath,1);
	  int_maskp(trpath(i,2),trpath(i,1)) = l+4*trpath(i,3);
	end
	
	imagesc(RtD*x,RtD*(y),int_maskp); set(gca,'ydir','normal');
	caxis([0,2*l])
	hold on
	PSIB = bstream(u*udim,zw*hdim,[y;ymax]*r0dim);
	img = PSIB(2:end,:)';
	contour(RtD*x,RtD*(y),img,40,'linewidth',.5);
	hold off
  end
  
  %% - Compute transport for path
  transport = compute_transport(trpath, u*udim, v*udim, ...
								[x;xmax]*cos(yv(2:m+1))'*r0dim, ...
								[y;ymax]*r0dim, zw*hdim);

  fprintf('transport through selected section: %f\n', transport);
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



