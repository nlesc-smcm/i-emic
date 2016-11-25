function [] = view_mask(mask_name)

  M = load([mask_name, '.mat']);

  fields = fieldnames(M);


  tmp_maskp = M.maskp;

  [m,n] = size(tmp_maskp);
  
  %% Display the transects already in there
  mx = max(max(abs(tmp_maskp)));
  basin_names = {};
  for i = 1:numel(fields)
 	if numel(fields{i}) == 2
	  fprintf('Found %s\n', fields{i});
	  pth = M.(fields{i});
	  for j = 1:size(pth,1);
		tmp_maskp(pth(j,2),pth(j,1)) = mx+1*pth(j,3);		
	  end
	end
	if numel(fields{i}) == 4
	  basin_names = [basin_names, fields{i}];
	end
  end

  figure(1)
  imagesc(tmp_maskp); set(gca,'ydir','normal');
  if numel(fields) > 1 
	for i = 1:numel(fields)
	  if numel(fields{i}) == 2
		pth = M.(fields{i});
		text(pth(end,1)+1, pth(end,2)+1, fields{i},'color','k','fontsize',10)
	  end
	end
  end

  cmap = [(7:.1:mx)'/mx,(7:.1:mx)'/mx,(7:.1:mx)'/mx];
  cmap = [cmap;1,0,0];
  colormap(cmap)
  
  xt = linspace(1,n,7);
  set(gca,'xtick',xt);
  set(gca,'xticklabels',linspace(0,360,numel(xt)))
  xlabel('Longitude')
  
  yt = linspace(1,m,9);
  set(gca,'ytick',yt);
  set(gca,'yticklabels',linspace(-80,80,numel(yt)))
  ylabel('Latitude')

  exportfig('transects.eps',9,[14,8])

  figure(2)
  cmap = [linspace(7,mx,mx+1)'/mx,linspace(7,mx,mx+1)'/mx,linspace(7,mx,mx+1)'/mx];
  cmap = [cmap;1,0,0];
  cmap = [cmap;0,0,1];
	
  colormap(cmap)
  %tmp_maskp = M.maskp;


  mx = max(max(tmp_maskp))
  plain_mask = mx* ones(m,n);
  
  for i = 1:numel(basin_names)
	plain_mask = plain_mask + 2*(M.(basin_names{i}) ~= 0);	  
  end
  max(max(plain_mask))
  
  imagesc(tmp_maskp);
  set(gca,'ydir','normal');
  hold on
  contour(plain_mask,1,'linewidth',1.5);
  hold off
  if numel(fields) > 1 
	for i = 1:numel(fields)
	  if numel(fields{i}) == 2
		pth = M.(fields{i});
		text(pth(end,1)+1, pth(end,2)+1, fields{i},'color','k','fontsize',10)
	  end
	end
  end
    xt = linspace(1,n,7);
  set(gca,'xtick',xt);
  set(gca,'xticklabels',linspace(0,360,numel(xt)))
  xlabel('Longitude')
  
  yt = linspace(1,m,9);
  set(gca,'ytick',yt);
  set(gca,'yticklabels',linspace(-80,80,numel(yt)))
  ylabel('Latitude')
end
