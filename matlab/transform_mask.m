function [peri,pmask] = transform_mask(mask_name)

  M        =  load([mask_name, '.mat']);
  mask     =  ~M.mask;
  [m,n,l]  =  size(mask);
  fid      =  fopen([mask_name, '.mask'],'w');

  % set padded mask pmask  
  pmask = ones(m+2,n+2,l+2);
  pmask(2:m+1,2:n+1,2:l+1) = mask;

  % find and set periodic boundaries
  for k = 1:l
	peri = mask(:,1,k) + mask(:,n,k);
	peri = logical([0; squeeze(logical(peri == 0)); 0]);
	pmask(peri, 1,   k+1) = 3;
	pmask(peri, n+2, k+1) = 3;
  end
  
  % write padded mask to file
  fprintf('  writing to %s\n',[mask_name, '.mask']);
  
  for k = l+2:-1:1
	fprintf(fid, '%% %1d %1d %1d %1d\n', n, m, l, l+3-k);
	for j = m+2:-1:1
	  for i = 1:n+2
		fprintf(fid,'%1d',pmask(j,i,k));
	  end
	  fprintf(fid, '\n');
	end
  end
  fclose(fid);
end
