function [] = plot_mode(V,mode)

  fprintf(1,'-------------- Plot EOF #%d-------------\n',mode)
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

  sol = zeros(nun, n, m, l);
  [qz,dfzt,dfzw] = gridstretch(zw);
  %% - READ SOLUTION --------------
  dim = n*m*(l+la)*nun;
  sol(:) = V(:,mode);
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
  figure(1)
  contourf(RtD*([y;ymax+dy/2]-dy/2),zw*hdim',PSIG',14);
  colorbar
  title('MOC')
  xlabel('latitude')
  ylabel('depth (m)')
  exportfig(['mstream_eof',num2str(mode),'.eps'],10,[20,7])

  figure(2)
  Sp2 = squeeze(mean(S,1)); 
  contourf(RtD*linspace(ymin,ymax,m),linspace(min(zw),max(zw),l)*hdim,Sl',14);
  colorbar
  title('Isohalines')
  xlabel('latitude')
  ylabel('depth (m)')
  exportfig(['isohalines_eof',num2str(mode),'.eps'],10,[20,7])

  figure(3)
  contourf(RtD*linspace(ymin,ymax,m),linspace(min(zw),max(zw),l)*hdim,Tl',14);
  colorbar
  title('Isothermals')
  xlabel('latitude')
  ylabel('depth (m)')
  exportfig(['isothermals_eof',num2str(mode),'.eps'],10,[20,7])

end
