fprintf('Usage: \n')
fprintf('  Left mouse:    set coordinate\n')
fprintf('  Right mouse:   retry\n')
fprintf('  Middle mouse:  abort\n')

SHARED_DIR = getenv('SHARED_DIR');

mask_path = [SHARED_DIR,'/i-emic/data/mkmask/paleo/',mask_name];

M = load([mask_path, '.mat']);
maskp     = M.maskp;
tmp_maskp = maskp;
[m,n] = size(maskp);
trs_maskp = zeros(m,n);
int_maskp = zeros(m,n);
N = max([n,m]);

coordinates = [];
l = max(max(maskp));
while size(coordinates,1) < 2

  imagesc(tmp_maskp); set(gca,'ydir','normal');
  
  [x,y,button] = ginput(1);
  x = round(x);
  y = round(y);

  if (maskp(y,x) > 0)
	fprintf('Only select land points!\n');
	continue;
  end
  
  if (button == 1)
	tmp_maskp(y,x) = l+8;
	coordinates = [coordinates;x,y];
  elseif (button == 3)
	continue;
  elseif (button == 2)
	break;
  end
end

imagesc(tmp_maskp); set(gca,'ydir','normal');
pause(.2);
assert(norm(size(coordinates)-[2,2])==0);

% Obtain linear function through coordinates
M = [coordinates(:,1),[1;1]]; b = coordinates(:,2);
sol =  M \ b;
A = sol(1); B = sol(2);

x1 = coordinates(1,1);
x2 = coordinates(2,1);

y1 = coordinates(1,2);
y2 = coordinates(2,2);


x = linspace(x1,x2,100*N);
fun = @(X) A*X+B;
list = []; dd = [];
x0 = x1;
y0 = fun(x0);
add_incr  = false;
for i = 1:numel(x)
  xr = round(x(i));
  y  = fun(x(i));
  yr = round(y);
  if (trs_maskp(yr,xr) == 0)
	trs_maskp(yr,xr) = 1;
	
	dy = y - y0;
	dx = x(i) - x0;

	if add_incr	
	  dd = [dd;[dx,dy]];
	end
	
	if maskp(yr,xr) > 0
	  int_maskp(yr,xr) = 1;
	  list = [list;xr,yr];
	  add_incr = true;
	else
	  add_incr = false;
	end	

	y0 = y;
	x0 = x(i);
  end
	 
end

imagesc(maskp + 100*trs_maskp); set(gca,'ydir','normal'); hold on
plot([x1,x2],[y1,y2],'r'); hold on

for i = 1:n
	plot([i-.5,i-.5],[.5,m+.5],'k');
end
for j = 1:m
  plot([.5,n+.5],[j-.5,j-.5],'k');
end

hold off
pause(.2);



disp([list,dd])
