function [] = remove_transect(mask_name);
  %% remove_transect(mask_name);
  %%
  %% Displays the mask and the transects, opens a dialog to
  %% remove transects.
		 
  M = load([mask_name, '.mat']);


  maskp = M.maskp;

  int_maskp = maskp;

  %% Display the transects  
  while true		
	fields = fieldnames(M);
	tmp_maskp = maskp;

	mx = max(max(abs(maskp)));
	for i = 1:numel(fields)
	  if numel(fields{i}) == 2
		fprintf('Found %s\n', fields{i});
		pth = M.(fields{i});
		for j = 1:size(pth,1);
		  tmp_maskp(pth(j,2),pth(j,1)) = mx+4*pth(j,3);		
		end
	  end
	end

	figure(1)
	imagesc(tmp_maskp); set(gca,'ydir','normal');

	if numel(fields) > 1 
	  for i = 1:numel(fields)
		if numel(fields{i}) == 2
		  pth = M.(fields{i});
		  text(pth(end,1)+1, pth(end,2), fields{i},'color','w')
		end
	  end  
	  pause(0.1);

	  char = input('Do you want to remove a transect? y/n ', 's');
	  if char == 'n'
		break
	  end

	  RField = input('Remove field: ', 's');
	  M = rmfield(M,RField);
	end
  end
  
  fields = fieldnames(M);
  fprintf('saving %d fields to %s.mat:\n', numel(fields), mask_name);
  for i = 1:numel(fields)
	fprintf('     %s\n', fields{i});
  end

  onOctave = (exist ('OCTAVE_VERSION', 'builtin') > 0);
  if onOctave	  
	save([mask_name,'.mat'], '-struct', '-mat' , 'M');
  else
	save([mask_name,'.mat'], '-struct', 'M');
  end
end
