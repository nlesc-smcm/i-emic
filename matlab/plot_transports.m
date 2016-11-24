mode = 1;

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

fnamelabels = {'-1', ...,
			   '-2', ...,
			   '-3', ...,
			   '-4', ...,
			   '-5', ...,
			   '-6', ...,
			   '-7', ...,
			   '-8', ...,
			   '-9', ...,
			   '-10', ...
			  };

N = 10;

MASK_PATH = [SHARED_DIR,'/i-emic/data/mkmask/'];

% Define transects:
transects = {'DR','IN','PA','SA','TA','TE'}

M = numel(transects)

transports = zeros(M,N);

if mode == 1
  for i = 1:N
	fname = [MASK_PATH,original_masks{i},'.mat'];
	fprintf('\n------------\nloading %s\n',fname);
	Mstruct = load(fname);
	for j = 1:M
	  trpath = Mstruct.(transects{j});
	  fprintf('\ncomputing %s transport using %s, %s\n', ...
			  transects{j}, statefiles{i}, datafiles{i});
	  transports(j,i) = compute_transports(statefiles{i}, datafiles{i}, trpath);
	end
  end

  plot(transports','.--','linewidth',1.0,'markersize',15)

  legend(transects,'location','northwest')
  set(gca,'xtick',1:N)
  set(gca,'xticklabels',labels)
  grid on
  ylabel('transport (Sv)')
  xlim([.5,N+.5]);

  fprintf('\n\n');
  exportfig('transports65to20.eps',10,[16,12]);
end

% plot all solutions
if mode == 2
  for i = 1:N  
	plot_ocean(statefiles{i},datafiles{i},labels{i},labels{i})
  end
  
end

% create interpolated plots for a movie
if mode == 3
  counter = 0;
  for i = 1:N
	for j = 1:4
	  counter = counter + 1
	  fnamelabel = ['-',num2str(counter)];
	  plot_ocean(statefiles{i},datafiles{i},labels{i},fnamelabel)
	end

	if i < N
	  for k = .1:.1:.9
		counter = counter + 1;
		fnamelabel = ['-',num2str(counter)];
		plot_ocean_interp(k,statefiles{i},statefiles{i+1},...
						  datafiles{i},datafiles{i+1},labels{i},fnamelabel);
	  end
	end
  end
end

