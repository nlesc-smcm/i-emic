%This is for a global 2deg mask with 12 vertical layers

if ~exist('mask65Ma') || force
  [peri,mask65Ma] = ...
  transform_mask('Mask_65Ma_lon1.5-3-358.5_lat-79.5-3-79.5',true);
end
if ~exist('mask45Ma') || force
  [peri,mask45Ma] = ...
  transform_mask('Mask_45Ma_lon1.5-3-358.5_lat-79.5-3-79.5',true);
end
if ~exist('mask25Ma') || force
  [peri,mask25Ma] = ...
  transform_mask('Mask_25Ma_lon1.5-3-358.5_lat-79.5-3-79.5',true);
end

lvl = 1; % level
mask = mask45Ma;

figure(1);
imagesc(mask(:,:,1+lvl));
set(gca,'ydir','normal')

if ~exist('A') || force
  A = mmread('jacobian.ocean');
end

[N,M] = size(A);

n = 120; m = 54; l = 12; nun = 6;

assert(N == n*m*l*nun);

% (u,v,w,p,T,S)
iu  = 1:6:N;
iv  = 2:6:N;
iw  = 3:6:N;
ip  = 4:6:N;

iuv = sort([iu,iv]);

Auv  = A(iuv,iuv);
Guv  = A(iuv,ip);
Duv  = A(ip, iuv);

if ~exist('iAuv') || force
  iAuv = inverseblockdiagonal(Auv,2);
end

C = -Duv * iAuv * Guv;

%force = true;
if ~exist('r') || force
  r = findzerorows(Duv);
end

RWS = load('zerorows');
TRB = zeros(n,m);

ctr = 0;
for j = 1:m
  for i = 1:n
	if mask(j+1,i+1,2)
	  continue;
	end 
	ctr = ctr + 1;
	TRB(i,j) = RWS(ctr);
  end
end

force  = false;
uvgrid = zeros(1,n,m,l);

uvgrid(r) = 1;

figure(3)
imagesc(squeeze(uvgrid(1,:,:,13-lvl))')
set(gca,'ydir','normal')
condest(C)

%for i = 1:12
%figure(i)
%MSK = squeeze(mask(2:end-1,2:end-1,1+i))';
%UVG = squeeze(uvgrid(1,:,:,13-i));
%imagesc((MSK+10*TRB)')
%set(gca,'ydir','normal')
%end
