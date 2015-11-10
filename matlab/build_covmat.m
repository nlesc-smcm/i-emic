function [Q,M,V,D,sol] = build_covmat(range, Q, ev)
  fprintf(1,'----------------------------------------------\n')
  %% - DEFINE CONSTANTS - ----------------------------------------------
  udim  = 0.1;                 %[m/s]   Velocity scale
  r0dim = 6.4e6;               %[m]     Radius of Earth
  T0    = 15;                  %[deg C] Reference temperature
  S0    = 0;                  %[psu]   Reference salinity
  RtD   = 180/pi;              %[-]     Radians to degrees
  %% - READ MASK - -----------------------------------------------------
  [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');
  surfm      = landm(2:n+1,2:m+1,l+1);  %Only interior surface points
  landm_int  = landm(2:n+1,2:m+1,2:l+1);
  dx         = (xu(n+1)-xu(1))/n;
  dy         = (yv(m+1)-yv(1))/m;
  dz         = (zw(l+1)-zw(1))/l;

  [qz,dfzt,dfzw] = gridstretch(zw);
  %% - READ SOLUTION --------------
  dim = n*m*(l+la)*nun;
  obs = numel(range);

  if nargin < 2 || isempty(Q)
	% data matrix
	M = read_many_forts(range, dim);

	% remove checkerboard modes in pressure
	p_idx   = 4:6:dim; % pressure indices
	[s1,s2] = checkerboard_modes(n,m,l);

	scale1   = s1(:)'*s1(:);
	scale2   = s2(:)'*s2(:);
	for i = 1:obs
	  M(i,p_idx) = M(i, p_idx) - (M(i, p_idx)*s1(:)/scale1)*s1(:)' ...
	                 - (M(i, p_idx)*s2(:)/scale2)*s2(:)';
	end
	
	% subtract mean to get centred data matrix
	mn = mean(M);
	M = M - repmat(mn,obs,1);
	Q = (1/(obs-1))*M'*M;
  else
	M = [];
  end
  sol = zeros(nun, n, m, l);
  
  if nargin < 3
	ev = 1;
  end
  
  fprintf(1,'------------- Calculating eigenvector--------\n')
  [V,D,flag] = eigs(Q,ev+4,'lr');
  fprintf(' flag = %d\n', flag);
  sol(:) = V(:,ev);
  fprintf(1,'-------------- Plot eigenvector -------------\n')

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

  %% build longitudinal average over non-land cells
  Sl = zeros(m,l);
  Tl = zeros(m,l);
  for k = 1:l
    for j = 1:m
      count = 0;
      for i=1:n            
        if landm_int(i,j,k) == 0
          count = count + 1;
          Sl(j,k) = Sl(j,k) + S(i,j,k);
          Tl(j,k) = Tl(j,k) + T(i,j,k);
        end
      end
      Sl(j,k) = Sl(j,k) / count;
      Tl(j,k) = Tl(j,k) / count;
    end
  end

  %%
  figure(2)
  contourf(RtD*[y;ymax],zw*hdim',PSIG',15);
  colorbar
  title('Overturning Streamfunction')
  xlabel('Latitude')
  ylabel('z (m)')
  exportfig('mstream.eps')

  figure(3)
  Sp = S(:,:,l);
  temp = flipud(Sp');
  %contourf(RtD*x,RtD*y,S0+Sp',15);
  plot(RtD*y, mean(Sp,1),'r');
  title('Surface Salinity');
  xlim(RtD*[ymin,ymax])
  xlabel('Latitude');


  figure(4)
  Sp2 = squeeze(mean(S,1)); 
  contourf(RtD*yv(1:end-1),z*hdim,Sl',15);
  colorbar
  title('Isohalines')
  xlabel('Latitude')
  ylabel('z (m)')
  exportfig('isohalines.eps')
  
  figure(5)
  Tp = T(:,:,l);
  plot(RtD*y, mean(Tp,1),'r');
  title('Surface Temperature');
  xlim(RtD*[ymin,ymax])
  xlabel('Latitude');

  figure(6)
  contourf(RtD*yv(1:end-1),z*hdim,Tl',15);
  colorbar
  title('Isothermals')
  xlabel('Latitude')
  ylabel('z (m)')
  exportfig('isohalines.eps')
  
  figure(7)
  Pp = p(:,:,l);
  plot(RtD*y, mean(Pp,1),'r');
  title('Surface Pressure');
  xlim(RtD*[ymin,ymax])
  xlabel('Latitude');

end
