function [] = plot_mask(landmask,x,y,varargin)

  RtD   = 180/pi; % Radians to degrees
  n     = numel(x);
  m     = numel(y);

  assert(size(landmask,1) == n, 'landmask not cropped')
  assert(size(landmask,2) == m, 'landmask not cropped')

  mask = landmask;
  if (size(landmask,3) > 1)
      mask = sum(landmask,3);
  end
  
  mask = mask / max(max(abs(mask)));
  mask = mask.^3;
  mask(mask<0) = 0;
  mask(mask>1) = 1;

  thres = .1;
  x     = RtD*x;
  y     = RtD*y;
  dx    = x(2)-x(1);
  dy    = y(2)-y(1);

  X      = zeros(5,n*m);
  Y      = zeros(5,n*m);
  CL     = zeros(n*m, 1, 3);
  CLEdge = zeros(5*n*m,3);

  col_ctr = 0;
  for i = 1:n
      for j = 1:m
          val = mask(i,j);
          if val >= thres
              col_ctr = col_ctr + 1;
              X(:,col_ctr) = [x(i)-dx/2; x(i)+dx/2; x(i)+dx/2; x(i)- ...
                              dx/2; x(i)-dx/2];
              Y(:,col_ctr) = [y(j)-dy/2; y(j)-dy/2; y(j)+dy/2; y(j)+ ...
                              dy/2; y(j)-dy/2];

              CL(col_ctr,1,:) = 1 - [val,val,val];
              CLEdge(5*(col_ctr-1)+1:5*col_ctr,:) = repmat(1-[val,val,val],5,1);
          end
      end
  end
  X  = X(:,1:col_ctr);
  Y  = Y(:,1:col_ctr);
  CL = CL(1:col_ctr,:,:);

  CLEdge = CLEdge(1:5*col_ctr,:,:);

  patch(X,Y,CL,'FaceVertexCData', CLEdge, 'EdgeColor', 'flat',varargin{:});

  %patch(X,Y,CL,'EdgeColor','none');

  %xlim([min(x),max(x)]);
  %ylim([min(y),max(y)]);

end
