function [coords, status] = getcoordinates(mask)
  %% getcoordinates(mask)
  %%
  %% use mouse input to generate coordinates
  %% for a straight zonal or meridional section
  
  reset_mask = mask;
  fprintf('\nUsage: \n')
  fprintf('  Left mouse:    set coordinate\n')
  fprintf('  Right mouse:   retry\n')
  fprintf('  Middle mouse:  quit\n')
  status = 0;
  coords = [];
  mx = max(max(mask));
  [m,n] = size(mask);
  x = -1; y = -1;
  while size(coords,1) < 2
	
	imagesc(mask); set(gca,'ydir','normal');
	hold on;
	plot([x-.5,x-.5],[.5,m+.5],'k');
	plot([x+.5,x+.5],[.5,m+.5],'k');
	plot([.5,n+.5],[y-.5,y-.5],'k');
	plot([.5,n+.5],[y+.5,y+.5],'k');
	hold off
	drawnow
		
	[x,y,button] = ginput(1);
	x = round(x);
	y = round(y);
	
	if (button == 1)
	  mask(y,x) = mx+8;
	  coords = [coords;x,y];
	elseif (button == 3)
	  mask = reset_mask;
	  coords = [];
	  continue;
	elseif (button == 2)
	  status = 1;
	  coords = [];
	  return;
	end

	if size(coords,1) == 2
	  d = coords(2,:) - coords(1,:);
	  if d(1) ~= 0 && d(2) ~= 0
		fprintf('We only allow straight meridional and zonal transects!\n')
		fprintf('Resetting...\n');
		mask = reset_mask;
		coords = [];
	  end
	end
  end  
  
  imagesc(mask); set(gca,'ydir','normal');
  pause(.2);
  
  assert(norm(size(coords)-[2,2])==0);
end
