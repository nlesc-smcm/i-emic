% Merge two masks into 9 intermediate masks

mask_name_1 = 'Mask_65Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8';
mask_name_2 = 'Mask_64Ma_lon1.5-3-358.5_lat-79.5-3-79.5_qz1.8';

periodic = true;

M1 = load([mask_name_1, '.mat']);
M2 = load([mask_name_2, '.mat']);

mask1 = M1.maskp;
mask2 = M2.maskp;

maskS = mask1;
for s = 0.1:0.1:0.9
  maskO = maskS;
  maskS = round(s*mask2 + (1-s)*mask1);


  figure(2)
  imagesc(maskS)
  set(gca,'ydir','normal')
  title(s)
  drawnow
  
  mdiff = max(max(abs(maskO - maskS)));
  fprintf('max dif = %d\n', mdiff)
  if (mdiff > 0)
	 fname = sprintf('%s/s%1.1f%s',mask_name_1,s,mask_name2);
	 maskp = maskS;
	 system(['mkdir -p ', mask_name_1]);
	 fprintf(' writing to %s\n', fname);
	 save([fname,'.mat'],'maskp')
	 transform_mask(fname, periodic);
  end
end
