shared_dir = [getenv('SHARED_DIR'), '/i-emic/data/mkmask/paleo/'];

mask_names = {'Mask_0Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8',
			  'Mask_5Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8',
			  'Mask_10Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8',
			  'Mask_15Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8',
			  'Mask_20Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8',
			  'Mask_25Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8',
			  'Mask_30Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8',
			  'Mask_35Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8',
			  'Mask_40Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8',
			  'Mask_45Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8',
			  'Mask_50Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8',
			  'Mask_55Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8',
			  'Mask_60Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8',
			  'Mask_65Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8' };

titles = {'0Ma',
		  '5Ma',
		  '10Ma',
		  '15Ma',
		  '20Ma',
		  '25Ma',
		  '30Ma',
		  '35Ma',
		  '40Ma',
		  '45Ma',
		  '50Ma',
		  '55Ma',
		  '60Ma',
		  '65Ma'};

N = numel(mask_names);
m = floor(sqrt(N));    % rows
n = floor(sqrt(N))+1;  % cols
for i = 1:min(m*n, N)	
  h = subplot(m, n, i);
  M = load([shared_dir,mask_names{i},'.mat']); 
  imagesc(M.maskp); set(gca,'ydir','normal'); colormap(gray);
  axis off;
  title(titles{i});

  posrow = get(h,'position');
  posrow(1) = mod(i-1,n)*0.17;
  posrow(2) = 0.7093-floor((i-1)/n)*0.25;
  set(h,'position',posrow)
  disp(posrow)
end

exportfig('allmasks.eps',10,[70,40]);
