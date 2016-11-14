function [path] = getpath(coords, mask)
  %%  getpath(coords, mask)
  %%  
  %%  given coordinates in a 2x2 matrix coords
  %%  return a path through the grid
  
  x1 = coords(1,1);
  y1 = coords(1,2);
  x2 = coords(2,1);
  y2 = coords(2,2);

  %% starting point
  x = x1;
  y = y1;
  
  %% difference with target
  dx = (x - x2) / max(abs(x - x2),1);
  dy = (y - y2) / max(abs(y - y2),1);

  %% initialize list that stores the path and orientation
  % orientation: or = 1 : u-direction
  %              or = 2 : v-direction

  or = 2;
  if dx == 0
	or = 1;
  end

  if 1 %(mask(y,x) > 0)
	path = [x,y,or];
  else
	path = [];
  end

  corner = false;
  while (dx ~= 0) || (dy ~= 0)

	%% traverse the transect, first in the x-direction
	if dx ~= 0
	  switch dx
		case 1
		  x = x - 1;
		case -1 
		  x = x + 1;
	  end
	  dx = (x - x2) / max(abs(x - x2),1);
	  or = 2;
	  if (dx == 0) && (dy ~= 0) 
		corner = true;
	  end

	  %% y-direction if x is aligned with x2
	elseif dy ~= 0
	  
	  switch dy
		case 1
		  y = y - 1;
		case -1 
		  y = y + 1;
	  end
	  dy = (y - y2) / max(abs(y - y2),1);
	  or = 1;
	  corner = false;
	end
	
	if ~corner && 1%(mask(y,x) > 0)
	  path = [path;x,y,or];
	end
  end
end
