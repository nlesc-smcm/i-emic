function [] = plot_mask(landmask,x,y)

  RtD   = 180/pi; % Radians to degrees  
  n     = numel(x);
  m     = numel(y);
  
  assert(size(landmask,1) == n, 'landmask not cropped')
  assert(size(landmask,2) == m, 'landmask not cropped')
  
  mask = sum(landmask,3);
  mask = mask / max(max(abs(mask)));
  mask = mask.^3;
  mask(mask<0) = 0;
  mask(mask>1) = 1;

  thres = .1;
  x     = RtD*x;
  y     = RtD*y;
  dx    = x(2)-x(1);
  dy    = y(2)-y(1);

  for i = 1:n
	for j = 1:m
	  val = mask(i,j);
	  if val > thres
		rectangle('Position',[x(i)-dx/2,y(j)-dy/2,dx,dy],'EdgeColor',1-[val,val,val],'FaceColor',1-[val,val,val]);
	  end
	end
  end

  xlim([min(x),max(x)]);
  ylim([min(y),max(y)]); 
  
end
