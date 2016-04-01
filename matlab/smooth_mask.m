function [mask, orig_mask] = smooth_mask(mask_name, method, tries)

  if nargin < 3
	tries = 4;
	if nargin < 2
	  method = 1;
	  if nargin < 1
		mask_name = 'natl_290_355_10_80';
	  end
	end
  end

  M = load([mask_name, '.mat']);
  mask = M.mask;
  
  [m,n,l] = size(mask);

  if method == 1
	n_hor_min = 2;
  elseif method == 2
	n_hor_min = 4;
  end
	
  n_ver_min = 1;

  sea  = 1;
  land = 0;

  figure(1)
  imagesc(mask(:,:,1));
  title('mask level 1');
  figure(1); figure(2);
  
  orig_mask = mask;
  
  for tr = 1:tries
	fix_ind = [];
	
	for k = 1:l
	  for j = 1:m
		for i = 1:n
		  if (mask(j,i,k) == sea)
			n_hor = count_hor_friends(i,j,k,n,m,l,mask,sea, method);
			if n_hor < n_hor_min
			  fix_ind = [fix_ind; i,j,k,land];
			end
		  else
			n_hor = count_hor_friends(i,j,k,n,m,l,mask,land,3);
			if n_hor < 3
			  fix_ind = [fix_ind; i,j,k,sea];
			end			
		  end
		end
	  end
	end

	for t = 1:size(fix_ind,1);
		
	  fi = fix_ind(t,1);
	  fj = fix_ind(t,2);
	  fk = fix_ind(t,3);
	  
	  if fix_ind(t,4) == land
		mask(fj, fi, fk:l) = land;
	  elseif fix_ind(t,4) == sea
		mask(fj, fi, fk) = sea;
	  end
	end

  end
  
  %  figure(2)
  %  imagesc(mask(:,:,1));
  %  title('mask level 1');
  %
  %  figure(3)
  %  imagesc(mask(:,:,5));
  %  title('mask level 5');
  %
  %  figure(4)
  %  imagesc(mask(:,:,10));
  %  title('mask level 10');

  save([mask_name,'.mat'], 'mask');

end

function [count] = count_hor_friends(i,j,k,n,m,l,mask,typ, method)
  nb = hor_neighbours(i,j,k,n,m,l, method);
  mt = mask(nb);
  mt = mt(mt == typ);
  count = numel(mt);
end

function [nb] = hor_neighbours(i,j,k,n,m,l, method)
  nb = zeros(8,1);

  nb(1) = find_mask_row(i+1,j,k,n,m,l);
  nb(2) = find_mask_row(i,j+1,k,n,m,l);
  nb(3) = find_mask_row(i,j-1,k,n,m,l);
  nb(4) = find_mask_row(i-1,j,k,n,m,l);

  if method == 2
	if nb(1) || nb(2)
	  nb(5) = find_mask_row(i+1,j+1,k,n,m,l);
	end
	if nb(2) || nb(4)
	  nb(6) = find_mask_row(i-1,j+1,k,n,m,l);
	end
	if nb(1) || nb(3)
	  nb(7) = find_mask_row(i+1,j-1,k,n,m,l);
	end
	if nb(3) || nb(4)
	  nb(8) = find_mask_row(i-1,j-1,k,n,m,l);
	end
  end

  if method == 3
	nb(5) = find_mask_row(i+1,j+1,k,n,m,l);
	nb(6) = find_mask_row(i-1,j+1,k,n,m,l);
	nb(7) = find_mask_row(i+1,j-1,k,n,m,l);
	nb(8) = find_mask_row(i-1,j-1,k,n,m,l);
  end
  
  nb = nb(nb>0);
end

function [row] = find_mask_row(i,j,k,n,m,l);
  if i < 1 || i > n || ...
	 j < 1 || j > m || ...
	 k < 1 || k > l
	row = 0;
  else
	row = j + (i-1)*m + (k-1)*m*n;
  end  
end
