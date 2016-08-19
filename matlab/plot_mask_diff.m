mask_name_1 = 'paleo/Mask_65Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8/coarsened';
mask_name_2 = 'paleo/Mask_64Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8/coarsened';

M1 = load([mask_name_1, '.mat']);
M2 = load([mask_name_2, '.mat']);

mask1 = M1.maskp;
mask2 = M2.maskp;

diff = mask1 - mask2;
figure(1); imagesc(mask1); set(gca,'ydir','normal');
figure(2); imagesc(mask2); set(gca,'ydir','normal');
figure(3); imagesc(diff);  set(gca,'ydir','normal');
