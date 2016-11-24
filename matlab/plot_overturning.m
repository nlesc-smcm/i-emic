create_data = true;

SHARED_DIR = getenv('SHARED_DIR');

statefiles = {'state_topo_0',...,
			  'state_topo_1',...,
			  'state_topo_2',...,
			  'state_topo_3',...,
			  'state_topo_4',...,
			  'state_topo_5',...,
			  'state_topo_6',...,
			  'state_topo_7',...,
			  'state_topo_8',...,
			  'state_topo_9',...
			 };

datafiles = {'mask_0.mask',...,
			 'mask_1.mask',...,
			 'mask_2.mask',...,
			 'mask_3.mask',...,
			 'mask_4.mask',...,
			 'mask_5.mask',...,
			 'mask_6.mask',...,
			 'mask_7.mask',...,
			 'mask_8.mask',...,
			 'mask_9.mask',...
			};

original_masks = {'paleo2/Mask_65Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8.dmp',...,
				  'paleo2/Mask_60Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8.dmp',...,
				  'paleo2/Mask_55Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8.dmp',...,
				  'paleo2/Mask_50Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8.dmp',...,
				  'paleo2/Mask_45Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8.dmp',...,
				  'paleo2/Mask_40Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8.dmp',...,
				  'paleo2/Mask_35Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8.dmp',...,
				  'paleo2/Mask_30Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8.dmp',...,
				  'paleo2/Mask_25Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8.dmp',...,
				  'paleo2/Mask_20Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8.dmp',...
				 };

labels = {'65Ma', ...,
		  '60Ma', ...,
		  '55Ma', ...,
		  '50Ma', ...,
		  '45Ma', ...,
		  '40Ma', ...,
		  '35Ma', ...,
		  '30Ma', ...,
		  '25Ma', ...,
		  '20Ma', ...
		 };

N = 10;

MASK_PATH = [SHARED_DIR,'/i-emic/data/mkmask/'];

								% Define basins:
basins = {'NOPA','SOPA','INOC','NOAT','SOAT'};

M = numel(basins);

[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = ...
readfort44(datafiles{1});

%% - DEFINE CONSTANTS - ----------------------------------------------
udim  = 0.1;                 %[m/s]   Velocity scale
r0dim = 6.4e6;               %[m]     Radius of Earth
T0    = 15;                  %[deg C] Reference temperature
S0    = 35;                  %[psu]   Reference salinity
RtD   = 180/pi;              %[-]     Radians to degrees

if create_data

  MOCdata = zeros(M,N,6);
  HOVdata = zeros(m+1,N,M);
  
  for i = 1:N
	fname = [MASK_PATH,original_masks{i},'.mat'];
	fprintf('\n------------\nloading %s\n',fname);
	Mstruct = load(fname);

	[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = ...
	readfort44(datafiles{i});
	dx         = (xu(n+1)-xu(1))/n;
	dy         = (yv(m+1)-yv(1))/m;
	dz         = (zw(l+1)-zw(1))/l;

	%% get the range in which we are interested
	%% in the z-direction everything below 1000m
	%% in the y-direction at latitudes 40N and 40S
	minDepth = 1000;
	latN   =  42;
	latS   = -38;
	
	zrange = (zw*hdim' < -minDepth);
	yrdim  = RtD*([y;ymax+dy/2] - dy/2);
	yrange = (abs(yrdim - latN) == min(abs(yrdim - latN))) + ...
			 (abs(yrdim - latS) == min(abs(yrdim - latS)));

	zrange = logical(zrange);
	yrange = logical(yrange);

	[lab icp par xl xlp det sig sol solup soleig] = ...
	readfort3(0, statefiles{i});

	%% - EXTRACT SOLUTION COMPONENTS - -----------------------------------
	[~,v,~,~,~,~] = extractsol(sol);
	
								% figure(1)
	PSIG = mstream(v*udim,[x;xmax]*cos(yv(2:m+1))'*r0dim,zw*hdim);
	PSIG = [zeros(m+1,1) PSIG];
			   % contourf(RtD*([y;ymax+dy/2]-dy/2),zw*hdim',PSIG',30);
			   % colorbar  

	for j = 1:M
	  v_restr   = v;
	  basinmask = Mstruct.(basins{j});
	  fprintf('  Basin: %s\n   ', basins{j});
	  basinmask = repmat(basinmask',[1 1 l]);
 	  v_restr(basinmask == 0) = 0;

								% figure(2);
								% imagesc(sum(v_restr,3)');
								% set(gca,'ydir','normal');
	  
								%figure(3)
	  PSIG = mstream(v_restr*udim,[x;xmax]*cos(yv(2:m+1))'*r0dim,zw*hdim);
	  PSIG = [zeros(m+1,1) PSIG];
	  PSIG(:,~zrange) = 0;
	  
								%PSIGplot = PSIG;
								%PSIGplot(:,~zrange) = NaN;
								%PSIG
	  
				%contourf(RtD*([y;ymax+dy/2]-dy/2),zw*hdim',PSIG',30);
				%colorbar

								% create Hovmoller data
	  HOVdata(:,i,j) = min(PSIG,[],2);
	  
	  psiMaxGlb = max(max(PSIG(:,zrange)));
	  psiMinGlb = min(min(PSIG(:,zrange)));

 	  if abs(psiMaxGlb) > abs(psiMinGlb)
		psiExtGlb = psiMaxGlb;
	  else
		psiExtGlb = psiMinGlb;
	  end	

	  psiMaxLat = max(max(PSIG(yrange,zrange)));
	  psiMinLat = min(min(PSIG(yrange,zrange)));

	  if abs(psiMaxLat) > abs(psiMinLat)
		psiExtLat = psiMaxLat;
	  else
		psiExtLat = psiMinLat;
	  end	

	  MOCdata(j,i,:) = [psiMinGlb,psiMaxGlb,...
						psiMinLat,psiMaxLat,...
						psiExtGlb,psiExtLat];
	end						   		
  end
end

figure(1)
plot(MOCdata(:,:,5)','.--','linewidth',1.0,'markersize',15)
legend(basins,'location','northwest')
set(gca,'xtick',1:N)
set(gca,'xticklabel',labels)
grid on
ylabel('MOC (Sv)')
title('Extremum in basin')
exportfig('overturningBasins.eps',10,[16,10]);

figure(2)
plot(MOCdata(:,:,6)','.--','linewidth',1.0,'markersize',15)
legend(basins,'location','northwest')
set(gca,'xtick',1:N)
set(gca,'xticklabel',labels)
grid on
ylabel('MOC (Sv)')
title('Extremum at latitude')
exportfig('overturningLatitudes.eps',10,[16,10]);

figure(3)
imagesc(1:N, RtD*([y;ymax+dy/2]-dy/2), HOVdata(:,:,2));
set(gca,'ydir','normal');
set(gca,'xtick',1:N)
set(gca,'xticklabel',labels)
colorbar

