% Merge two masks into X intermediate masks

shared_dir = [getenv('SHARED_DIR'), '/i-emic/data/mkmask/paleo/'];

mask_name_1 = 'Mask_55Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8';
mask_name_2 = 'Mask_50Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8';

periodic   = true;
savemerged = true;

M1 = load([shared_dir, mask_name_1, '.mat']);
M2 = load([shared_dir, mask_name_2, '.mat']);

mask1 = M1.maskp;
mask2 = M2.maskp;

maskS = mask1;
srange = 0.1:0.1:0.9;
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
  
  if (mdiff(end) > 0)
	if (savemerged)
	  fname  = sprintf('%s%s/s%1.1f%s',shared_dir,mask_name_1,s,mask_name_2);
	  maskp  = maskS;
	  newdir = sprintf('%s%s', shared_dir, mask_name_1);
	  fprintf('create new dir: %s\n', newdir);
	  system(['mkdir -p ', newdir])
	  fprintf(' writing to %s\n', fname);
	  save([fname,'.mat'],'maskp')
	  transform_mask(fname, periodic);
	end
  end
input('press ENTER to continue');
end
plot(srange, ndiff)

figure(2)
imagesc(mask2)
set(gca,'ydir','normal')
title(s)
drawnow
