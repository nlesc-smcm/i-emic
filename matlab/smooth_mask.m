function [mask] = smooth_mask(mask_name)
  if nargin < 1
	mask_name = 'natl_290_355_10_80';
  end

  M = load([mask_name, '.mat']);
  mask = M.mask;
  
  [m,n,l] = size(mask);
  
  n_hor_min = 2;
  n_ver_min = 1;

  sea  = 1;
  land = 0;

  figure(1)
  imagesc(mask(:,:,1));
  title('mask level 1');

  i = 15; j = 10; k = 1;
  mask(j,i,k)
  n_hor = count_hor_friends(i,j,k,n,m,l,mask,sea);

  fix_ind = [];

  tries = 8;
  
  for t = 1:tries
  
  for k = 1:l
	for j = 1:m
	  for i = 1:n
		if (mask(j,i,k) == sea)
		  n_hor = count_hor_friends(i,j,k,n,m,l,mask,sea);
		  if n_hor < n_hor_min
			fix_ind = [fix_ind; i,j,k];
		  end
		end
	  end
	end
  end

  for t = 1:size(fix_ind,1);
	  
	fi = fix_ind(t,1);
	fj = fix_ind(t,2);
	fk = fix_ind(t,3);
	
	mask(fj, fi, fk:l) = land;

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

function [count] = count_hor_friends(i,j,k,n,m,l,mask,sea)
  nb = hor_neighbours(i,j,k,n,m,l);
  mt = mask(nb);
  mt = mt(mt == sea);
  count = numel(mt);
end

function [nb] = hor_neighbours(i,j,k,n,m,l)
		 nb = zeros(4,1);
		 
		 nb(1) = find_mask_row(i+1,j,k,n,m,l);
		 nb(2) = find_mask_row(i,j+1,k,n,m,l);
		 nb(3) = find_mask_row(i,j-1,k,n,m,l);
		 nb(4) = find_mask_row(i-1,j,k,n,m,l);
		 
		 nb = nb(nb > 0);
end

function [row] = find_mask_row(i,j,k,n,m,l);
  if i < 1 || i > n || ...
	 j < 1 || j > m || ...
	 k < 1 || k > l
	row = -1;
  else
	row = j + (i-1)*m + (k-1)*m*n;
  end  
end
