% Merge two masks into X intermediate masks

mask_name = 'Mask_65Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8';

periodic   = true;
savemerged = false;

M = load([mask_name, '.mat']);

maskp = M.maskp;

figure(1); imagesc(maskp); set(gca,'ydir','normal');
title('undampened')
maskp(maskp==1) = 2;
figure(2); imagesc(maskp); set(gca,'ydir','normal');
title('dampened')

if (savemerged)
  fname = sprintf('%s.dmp', mask_name);
  fprintf(' writing to %s\n', fname);
  save([fname,'.mat'],'maskp')
  transform_mask(fname, periodic);
end
