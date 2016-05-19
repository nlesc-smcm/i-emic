%This is for a global 2deg mask with 12 vertical layers

[peri,mask] = transform_mask('Mask_65Ma_lon1.5-3-358.5_lat-79.5-3-79.5',true);

figure(1);
imagesc(mask(:,:,3));
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

figure(2)
spy(iAuv(1:1000,1:1000))
force = false;

C = -Duv*iAuv*Guv;
condest(C)

%r = findzerorows(diag(C));

uvgrid = zeros(1,n,m,l);

uvgrid(r) = 1;

imagesc(squeeze(uvgrid(1,:,:,11))')
set(gca,'ydir','normal')
