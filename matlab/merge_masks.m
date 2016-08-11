% Merge two masks into X intermediate masks

mask_name_1 = 'Mask_40Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8';
mask_name_2 = 'Mask_30Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8';

periodic   = true;
savemerged = false;

M1 = load([mask_name_1, '.mat']);
M2 = load([mask_name_2, '.mat']);

mask1 = M1.maskp;
mask2 = M2.maskp;

maskS = mask1;
srange = 0:0.2:1;
mdiff = [];
ndiff = [];
for s = srange
  maskO = maskS;
  maskS = round((s*mask2) + ((1-s)*mask1));


  figure(2)
  imagesc(maskS)
  set(gca,'ydir','normal')
  title(s)

  figure(3)
  imagesc(maskS-maskO)
  set(gca,'ydir','normal')
  title(s)  
  
  mdiff = [mdiff, max(max(abs(maskO - maskS)))];
  ndiff = [ndiff, nnz(maskO - maskS)];
  fprintf('max dif = %d\n', mdiff(end))
  fprintf('num dif = %d\n', ndiff(end))
  
  pause(0.2)
  if (mdiff > 0)
	if (savemerged)
	  fname = sprintf('%s/s%1.1f%s',mask_name_1,s,mask_name_2);
	  maskp = maskS;
	  system(['mkdir -p ', mask_name_1]);
	  fprintf(' writing to %s\n', fname);
	  save([fname,'.mat'],'maskp')
	  transform_mask(fname, periodic);
	end
  end
end
plot(srange, ndiff)

figure(2)
imagesc(mask2)
set(gca,'ydir','normal')
title(s)
drawnow
