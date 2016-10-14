function [] = edit_mask(mask_name)

  M = load([mask_name, '.mat']);
  maskp    = M.maskp;
  dmp_name = [mask_name,'.dmp'];  
  
  fprintf('Usage: \n')
  fprintf('  Left mouse:    decrease depth\n')
  fprintf('  Right mouse:   increase depth\n')
  fprintf('  Middle mouse:  quit\n')

  new_maskp = isolate_basin(maskp);
  imagesc(new_maskp); set(gca,'ydir','normal'); 

  l = max(max(maskp));

  while true		
	while true
		  
	  [x,y,button] = ginput(1);
	  x = round(x);
	  y = round(y);
	  
	  if (button == 1)
		new_maskp(y,x) = min(new_maskp(y,x)-1,l);
	  elseif (button == 3)
		new_maskp(y,x) = max(new_maskp(y,x)+1,0);
	  elseif (button == 2)
		break;
	  end
	  
	  new_maskp = isolate_basin(new_maskp);
	  imagesc(new_maskp); set(gca,'ydir','normal');	
	end

	maskdiff = new_maskp - maskp;
	imagesc(maskdiff); set(gca,'ydir','normal');

	title('difference')

	char = input('Satisfied? y/n ', 's');

	if (char == 'y')
	  maskp = new_maskp;
	  save([dmp_name, '.mat'], 'maskp');
	  transform_mask(dmp_name, true);
	  break;
	else
	  imagesc(new_maskp); set(gca,'ydir','normal');	
	end
	
  end
  
end

