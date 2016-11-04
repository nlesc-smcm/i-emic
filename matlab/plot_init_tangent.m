function [] = plot_init_tangent(fname, maskname)

  if nargin < 1
	fname = 'initial_tangent_backup.ocean';
  end
  if nargin < 2
	maskname = 'mask_1.mask';
  end

  [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44(maskname);

  RtD   = 180/pi;   %[-]     Radians to degrees

  res  = load(fname);

  tol = 1e-7;
  res(res>tol)  =  1;
  res(res<-tol)  = -1;

  sol = zeros(nun,n,m,l+la);
  idx = 1;
  for k = 1:l+la
	for j = 1:m
	  for i = 1:n
		for XX = 1:nun
		  sol(XX,i,j,k) = res(idx);
		  idx = idx + 1;
		end
	  end
	end
  end

  [u,v,w,p,T,S] = extractsol(sol);

  level = l;

  land = landm(2:end-1,2:end-1,2:l+1);
  
  cmap = [0,0,.85;
		  1,1,1;
	      .85,0,0;
		  .85,.85,.85];

  figure(1);
  img = 2*land(:,:,level)' + u(:,:,level)';
  imagesc(RtD*x,RtD*(y), img);
  set(gca,'ydir','normal');
  
  colormap(cmap)
  xlabel('Longitude')
  ylabel('Latitude');
  %colorbar
  title('u');
  exportfig('tangent_u.eps',10,[15,10]);
  
  figure(2);
  img = 2*land(:,:,level)' + v(:,:,level)';
  imagesc(RtD*x,RtD*(y), img);
  %colorbar
  title('v')
  set(gca,'ydir','normal');

  colormap(cmap);
  xlabel('Longitude');
  ylabel('Latitude');

  exportfig('tangent_v.eps',10,[15,10]);
end
