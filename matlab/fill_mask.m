function [] = fill_mask(mask_name)
		 
  dmp_name = [mask_name,'.dmp'];
  convert_to_txt(mask_name, dmp_name);
  
  fprintf('Edit the file in your favorite editor, then press ENTER...\n');
  pause(1);
  input('press ENTER');

  convert_to_bin(dmp_name);

  remove_inland_seas(dmp_name);

  transform_mask(dmp_name, true);

end

