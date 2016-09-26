shared_dir = [getenv('SHARED_DIR'), '/i-emic/data/mkmask/paleo/'];

mask_name_1 = 'Mask_20Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8';
mask_name_2 = 'Mask_15Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8';

M1 = load([shared_dir, mask_name_1, '.mat']);
M2 = load([shared_dir, mask_name_2, '.mat']);

mask1 = M1.maskp;
mask2 = M2.maskp;

diff = mask1 - mask2;

figure(8); imagesc(mask1);  title(mask_name_1,'interpreter', 'none'); set(gca,'ydir','normal'); grid on; colorbar
figure(9); imagesc(mask2);  title(mask_name_2,'interpreter', 'none'); set(gca,'ydir','normal'); grid on; colorbar
figure(10); imagesc(diff);   title('diff');      set(gca,'ydir','normal'); grid on; colorbar
