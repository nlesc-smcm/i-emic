function [] = plot_failed_residual(fname, maskname, level)

  if nargin < 1
	fname = 'residual.ocean';
  end
  if nargin < 2
	maskname = 'fort.44';
  end
  if nargin < 3
	level = l;
  end

  [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44(maskname);

  RtD   = 180/pi;              %[-]     Radians to degrees
  

  res  = load(fname);

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


  land = landm(2:end-1,2:end-1,2:l+1); 

  figure(1);
  mx = max(max(u(:,:,level)'));
  img = -mx*2*land(:,:,level)' + u(:,:,level)';
  imagesc(RtD*x,RtD*(y), img);
  set(gca,'ydir','normal');

  colormap(parula)
  xlabel('Longitude')
  ylabel('Latitude');
  colorbar
  title('u')  

  figure(2);
  mx = max(max(v(:,:,level)'));
  img = -mx*2*land(:,:,level)' + v(:,:,level)';
  imagesc(RtD*x,RtD*(y), img);
  colorbar
  title('v')
  set(gca,'ydir','normal');
  colormap(parula)
  xlabel('Longitude')
  ylabel('Latitude');

  figure(3);
  mx = max(max(w(:,:,level)'));
  img = -mx*2*land(:,:,level)' + w(:,:,level)';
  imagesc(RtD*x,RtD*(y), img);
  colorbar
  title('w')
  set(gca,'ydir','normal');
  colormap(parula)
  xlabel('Longitude')
  ylabel('Latitude');

  figure(4);
  mx  = max(max(p(:,:,level)'));
  img = -mx*2*land(:,:,level)' + p(:,:,level)';
  imagesc(RtD*x,RtD*(y), img);
  colorbar
  title('p')
  set(gca,'ydir','normal');
  colormap(parula)
  xlabel('Longitude')
  ylabel('Latitude');

  figure(5);
  mx  = max(max(T(:,:,level)'));
  img = -mx*2*land(:,:,level)' + T(:,:,level)';
  imagesc(RtD*x,RtD*(y), img);
  colorbar
  title('T')
  set(gca,'ydir','normal');
  colormap(parula)
  xlabel('Longitude')
  ylabel('Latitude');

  figure(6);
  mx  =  max(max(S(:,:,level)'));
  img = -mx*2*land(:,:,level)' + S(:,:,level)';
  imagesc(RtD*x,RtD*(y), img);
  colorbar
  title('S')
  set(gca,'ydir','normal');
  colormap(parula)
  xlabel('Longitude')
  ylabel('Latitude');

  fprintf('|u| = %e\n',norm(u(:)));
  fprintf('|v| = %e\n',norm(v(:)));
  fprintf('|w| = %e\n',norm(w(:)));
  fprintf('|p| = %e\n',norm(p(:)));
  fprintf('|T| = %e\n',norm(T(:)));
  fprintf('|S| = %e\n',norm(S(:)));
  fprintf('---------\n');
  fprintf('|X| = %e\n',norm(res(:)));
end
