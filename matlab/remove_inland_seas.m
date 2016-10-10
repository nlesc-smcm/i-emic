function [] = remove_inland_seas(mask_name)
  global cntr masktmp
  
  M = load([mask_name, '.mat']);
  maskp = M.maskp;

  [m,n] = size(maskp);
  masktmp = maskp;

  level = 1; % division level default=1
  
  masktmp(masktmp < level) = 0;
  masktmp(masktmp >= level) = 1;
  
  % start somewhere in the pacific
  % j = floor(m/2);
  % i = floor(1.7*n/3);

  % color everything that is reachable under right angles
  cntr = 0;
  
  %  cntr = 1;
  %  for dr = 1:4
  % 	color(dr, i, j, m, n);
  %  end


  maskperiod = [masktmp(:,round(n/2):end),masktmp,masktmp(:,1:round(n/2))];
  id  = numel(round(n/2):n)+1:numel(round(n/2):n)+n;  

  CC = bwconncomp(maskperiod,4);
  
  id1 = CC.PixelIdxList{1}; % this is tricky
  
  maskperiod(id1) = 999;  
  maskperiod(maskperiod~=999) = 0;
  maskp(maskperiod(:,id)~=999) = 0;

  fprintf(' removing inland seas\n');
  fprintf(' saving to %s.mat\n', mask_name);
  save([mask_name,'.mat'], 'maskp');
  
  imagesc(maskp);
  set(gca,'ydir','normal');
  drawnow
  
end


function [] = color(direction, i, j, m, n)
  global cntr masktmp


  direction
  cntr = cntr + 1;
  switch direction
	case 1
	  i = mod(i+1, n);
	case 3
	  i = mod(i-1, n);
	case 2
	  j = max(min(j + 1, m),1);
	case 4
	  j = max(min(j - 1, m),1);
  end

  if i == 0
	i = n
  end
  
  if masktmp(j,i) == 0 || masktmp(j,i) == 3
	return;
  end
  
  masktmp(j,i) = 3;

  if mod(cntr,50) == 0
	imagesc(masktmp);
	set(gca,'ydir','normal');
	drawnow
  end
  
  for dr = 1:4
	color(dr,i,j,m,n)
  end
  
end
