function [] = plot_ocean(solfile, maskfile, title_add, fname_add)
%---------------------------------------------------------------------
% PLOTTHCM - Mother script for plotting THCM output
%  usage: plot_ocean(solfile, datafile, title_add, fname_add)
%
%  Father is M. den Toom, who conceived it 06-11-08     
%  Modified by Erik, 2015/2016/2017 -> t.e.mulder@uu.nl
%---------------------------------------------------------------------

  title_additional = '';
  fname_additional = '';
  specify_mask = true; 

  if nargin < 1
	solfile = 'fort.3';
  end
  if nargin < 2
	maskfile = 'fort.44';
	specify_mask = false;
  end
  if nargin >= 3
	title_additional = title_add;
	fname_additional = fname_add;
  end

  fprintf(1,'----------------------------------------------\n')

  %% - DEFINE CONSTANTS - ----------------------------------------------

  udim  = 0.1;                 %[m/s]    Velocity scale
  r0dim = 6.4e6;               %[m]      Radius of Earth
  T0    = 15;                  %[deg C]  Reference temperature
  S0    = 35;                  %[psu]    Reference salinity
  RtD   = 180/pi;              %[-]      Radians to degrees

  %% - READ MASK - -----------------------------------------------------

  [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = ...
  readfort44(maskfile);

  surfm      = landm(2:n+1,2:m+1,l+1);  %Only interior surface points
  landm_int  = landm(2:n+1,2:m+1,2:l+1);
  dx         = (xu(n+1)-xu(1))/n;
  dy         = (yv(m+1)-yv(1))/m;
  dz         = (zw(l+1)-zw(1))/l;

 % - Create surface landmask image
  srf = [];
  greyness = 1;
  summask = sum(landm_int,3);
  summask = summask / max(max(abs(summask)));
  summask = summask.^3;
  
  srf(:,:,1) = (1-greyness*((summask')));
  srf(:,:,2) = (1-greyness*((summask')));
  srf(:,:,3) = (1-greyness*((summask')));
  srf(srf<0) = 0;
  srf(srf>1) = 1;

								% - Deduce grid stretching
  [qz,dfzt,dfzw] = gridstretch(zw);

  %% - READ SOLUTION - -------------------------------------------------
  [lab icp par xl xlp det sig sol solup soleig] = ...
  readfort3(la, solfile);

  %% - EXTRACT SOLUTION COMPONENTS - -----------------------------------
  [u,v,w,p,T,S] = extractsol(sol);

  %% - INTEGRATIONS - --------------------------------------------------
  % Barotropic streamfunction;
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
  Sl = zeros(m,l);

  for k = 1:l
	for j = 1:m
      count = 0;
      for i=1:n
		if landm_int(i,j,k) == 0
          count = count + 1;
          Tl(j,k) = Tl(j,k) + T(i,j,k);
		  Sl(j,k) = Sl(j,k) + S(i,j,k);
		else
		  S(i,j,k) = 0;
		  T(i,j,k) = 0;
		end
      end
      Tl(j,k) = Tl(j,k) / count;
	  Sl(j,k) = Sl(j,k) / count;
	end
  end

  % --- Create colormaps
  par = [0    0.4470    0.7410;  0.8500    0.3250    0.0980];
  neg  = par(1,:); 
  pos  = par(2,:); 
  N1   = 64;
  N2   = 64;
  mid  = [1,1,1];
  col1 = [linspace(neg(1),mid(1),N1)',linspace(neg(2),mid(2),N1)',linspace(neg(3),mid(3),N1)'];
  col2 = [linspace(mid(1),pos(1),N2)',linspace(mid(2),pos(2),N2)',linspace(mid(3),pos(3),N2)'];
  col_white  = [col1;col2];

  mid  = [1,1,1];
  col1 = [linspace(neg(1),mid(1),N1)',linspace(neg(2),mid(2),N1)',linspace(neg(3),mid(3),N1)'];
  col2 = [linspace(mid(1),pos(1),N2)',linspace(mid(2),pos(2),N2)',linspace(mid(3),pos(3),N2)'];
  col_black  = [col1;col2];



  %% - PLOT THE RESULTS - ----------------------------------------------
  figure(1)
  img  = PSIB(2:end,:)';
  imgp = img; imgp(imgp<0)=NaN;
  imgn = img; imgn(imgn>-0)=NaN;

  contourf(RtD*x,RtD*(y),imgp,15,'LineStyle','none'); hold on;
  contourf(RtD*x,RtD*(y),imgn,15,'LineStyle','none'); hold on;

  plot_mask(summask,x,y); hold on

  contour(RtD*x,RtD*(y),imgp,15,'k'); hold on;
  contour(RtD*x,RtD*(y),imgn,15,'k--'); hold off;
  
  colorbar
  colormap(col_black)
  %title(['Barotropic Streamfunction (Sv) ', title_additional]);
  xlabel('Longitude')
  ylabel('Latitude');
  %caxis([-40,40]);
  exportfig(['bstream',fname_additional,'.eps'],10,[19,11])

%%%
  figure(2)
  PSIGp = PSIG; PSIGp(PSIGp<0) = NaN;
  PSIGn = PSIG; PSIGn(PSIGn>0)  = NaN;
  contourf(RtD*([y;ymax+dy/2]-dy/2),zw*hdim',PSIGp',15); hold on
  contourf(RtD*([y;ymax+dy/2]-dy/2),zw*hdim',PSIGn',15,'--'); hold off
  colorbar
  cmin = min(min(PSIG(:,1:9)))
  cmax = max(max(PSIG(:,1:9)))
  %title(['MOC (Sv) ',title_additional])
  xlabel('latitude')
  ylabel('depth (m)')
  %caxis([-40,40])
  colormap(col_white)

  exportfig(['mstream',fname_additional,'.eps'],10,[19,10])
  

  %% -------------------------------------------------------
  figure(4)
  contourf(RtD*yv(1:end-1),z*hdim,Tl'+T0,15);
							  %imagesc(RtD*yv(1:end-1),z*hdim,Tp'+T0);
							  %set(gca,'ydir','normal');
							  %pcolor(RtD*yv(1:end-1),z*hdim,Tp'+T0);
  colorbar
  title('Temperature')
  xlabel('Latitude')
  ylabel('z (m)')
  colormap(col_white)
  exportfig('isothermals.eps',10,[20,7])
  return

  %% -------------------------------------------------------
  figure(3)
  Tsurf = T(:,:,l);
  minT = T0+min(min(Tsurf));
  maxT = T0+max(max(Tsurf));

  img  = T0 + Tsurf';
  contourf(RtD*x,RtD*(y),img,20,'Visible', 'off'); hold on;
  set(gca,'color',[0.65,0.65,0.65]);
					%image(RtD*x,RtD*(y),srf,'AlphaData',0.5); hold on
  contours = linspace(minT,maxT,20);
  imagesc(RtD*x,RtD*(y),img);
								%imagesc(RtD*x,RtD*(y),img)
								%set(gca,'ydir','normal');
  hold off

  colorbar
  title('Surface Temperature', 'interpreter', 'none');
  xlabel('Longitude');
  ylabel('Latitude');
  colormap(col_white)
  exportfig('sst.eps',10,[50,25])


  %%
  figure(6)
  contourf(RtD*yv(1:end-1),z*hdim,Sl'+S0,15);
							   %imagesc(Sp'+S0);
							   %set(gca,'ydir','normal')
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
  for j = 1:m
	for i = 1:n
	  if surfm(i,j) == 1
		Ssurf(i,j) = -999;
	  end
	end
  end

  img  = S0 + Ssurf';
  contourf(RtD*x,RtD*(y),img,20,'Visible', 'off'); hold on;
  set(gca,'color',[0.65,0.65,0.65]);
  image(RtD*x,RtD*(y),srf,'AlphaData',0.5); hold on
  contours = linspace(minS,maxS,40);
  contourf(RtD*x,RtD*(y),img,contours,'Visible', 'on','linewidth',1); hold off

								%imagesc(RtD*x,RtD*(y),img);
  set(gca,'ydir','normal');
  grid on
								%caxis([minS,maxS]);
  colorbar
  title('Surface Salinity', 'interpreter', 'none');
  xlabel('Longitude');
  ylabel('Latitude');


  exportfig('sss.eps',10,[50,25])

end
