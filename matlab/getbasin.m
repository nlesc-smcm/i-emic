function [mask_out] = getbasin(coords, mask_in)

  [m,n] = size(mask_in);
  
  %% obtain ordered coordinates defining the rectangle
  %% such that x1 < x2, y1 < y2  
  x1 = min(coords(1,1),coords(2,1));
  x2 = max(coords(1,1),coords(2,1));
  y1 = min(coords(1,2),coords(2,2));
  y2 = max(coords(1,2),coords(2,2));

  periodic = false;
  if (n+x1-x2) < (x2-x1)
	char = input('Does the basin cross a periodic boundary? y/n ', 's');
	if char == 'y'
	  periodic = true;
	end
  end	
  
  mask_out = mask_in;

  %% exclude everything except the interior and the boundary
  if ~periodic 
	mask_out(1:y1-1,:)   = 0;
	mask_out(y2+1:end,:) = 0;
	mask_out(:,1:x1-1)   = 0;
	mask_out(:,x2+1:end) = 0;
  elseif periodic
	mask_out(1:y1-1,:)    = 0;
	mask_out(y2+1:end,:)  = 0;
	mask_out(:,x1+1:x2-1) = 0;
  end

  % this needs some work
  %mask_out = isolate_basin(mask_out);
end
