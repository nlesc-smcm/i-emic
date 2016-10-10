function [peri, pmask] = transform_mask(mask_name, periodic)
  % This function transforms Michiel's masks to THCM masks
  % Supply the name of the .mat file and whether you want periodic
  % boundaries. 
		 
  if nargin < 2
	periodic = true;
  end

  M        =   load([mask_name, '.mat']);
  mimport  =   M.maskp; % might be .mask or .maskp, not sure
  [~,~,L]  =   size(mimport);

  reducesize = false;
  if (reducesize && L == 1) % cruel reduction
	mimport = mimport(1:2:end,1:2:end);
  end
  
  mask     =  ~mimport;
  [m,n,l]  =  size(mask)

  % check whether we have a non-layered format
  if (l == 1)
	fprintf('  not in layered format: converting!\n');
	% convert to layered format
	l = max(max(abs(mimport)))
	mask = zeros(m,n,l);
	for i = 1:l
	  mask(:,:,i) = ~(mimport >= i);
	end
  end

  % Create file id 
  fid      =  fopen([mask_name, '.mask'],'w');
  
  % set padded mask pmask  
  pmask = ones(m+2,n+2,l+2);
  pmask(2:m+1,2:n+1,2:l+1) = mask;

  if (periodic)
	% find and set periodic boundaries
	for k = 1:l
	  peri = mask(:,1,k) + mask(:,n,k);
	  peri = logical([0; squeeze(logical(peri == 0)); 0]);
	  pmask(peri, 1,   k+1) = 3;
	  pmask(peri, n+2, k+1) = 3;
	end
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
