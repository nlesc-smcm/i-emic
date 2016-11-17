SHARED_DIR = getenv('SHARED_DIR');

statefiles = {'state_topo_0',...,
			  'state_topo_1',...,
			  'state_topo_2',...,
			  'state_topo_3',...
			 };

datafiles = {'mask_0.mask',...,
			 'mask_1.mask',...,
			 'mask_2.mask',...,
			 'mask_3.mask',...
			};

original_masks = {'paleo2/Mask_65Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8.dmp',...,
				  'paleo2/Mask_60Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8.dmp',...,
				  'paleo2/Mask_55Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8.dmp',...,
				  'paleo2/Mask_50Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8.dmp',...
				 };

xlabels = {'65Ma', ...,
		   '60Ma', ...,
		   '55Ma', ...,
		   '50Ma' ...
		  }

N = numel(xlabels);

MASK_PATH = [SHARED_DIR,'/i-emic/data/mkmask/'];

% Define transects:
transects = {'DR','IN','PA','SA','TA','TE'}

M = numel(transects)

transports = zeros(M,N);

for i = 1:N
  fname = [MASK_PATH,original_masks{i},'.mat'];
  fprintf('loading %s\n',fname);
  Mstruct = load(fname);
  for j = 1:M
	trpath = Mstruct.(transects{j});
	fprintf('computing %s transport using %s, %s\n', ...
			transects{j}, statefiles{i}, datafiles{i});
	transports(j,i) = compute_transports(statefiles{i}, datafiles{i}, trpath);
  end
end

plot(transports','.--','linewidth',1.5,'markersize',20)

legend(transects,'location','northwest')
set(gca,'xtick',1:N)
set(gca,'xticklabels',xlabels)
grid on
ylabel('transport (Sv)')
xlim([.5,N+.5]);
