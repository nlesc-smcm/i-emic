%Auv   = load('Auv');   Auv   = spconvert(Auv);
%BDAuv = load('BDAuv'); BDAuv = spconvert(BDAuv);
%CHAT  = load('CHAT');  CHAT  = spconvert(CHAT);
%TMP   = load('TMP');   TMP   = spconvert(TMP);
%CHAT2 = load('CHAT2'); CHAT2 = spconvert(CHAT2);
%

stop_here = true;
create_mask
stop_here = false;

[mask1, orig_mask] = smooth_mask(mask_name, 1, 1);
[mask2, ~]         = smooth_mask(mask_name, 1, 10);

[m,n,l] = size(orig_mask);

mask1 = abs(mask1 - 1);
mask2 = abs(mask2 - 1);
orig_mask = abs(orig_mask - 1);

figure(1);
for i = 1:l
  subplot(3,4,i)
  colormap jet
  imagesc(flipud(orig_mask(:,:,i) + 2*mask1(:,:,i) + 4*mask2(:,:,i)));
  title(i)
  
end
