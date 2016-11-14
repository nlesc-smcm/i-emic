function [] = append_transect(mask_name, remove_field);
  %% append_transect(mask_name);
  %%
  %% Append a transect to mask mask_name

  if nargin < 2
	remove_field = false;
  end
  
  M = load([mask_name, '.mat']);

  fields = fieldnames(M);

  maskp = M.maskp;

  tmp_maskp = maskp;
  int_maskp = maskp;

  %% Display the transects already in there
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

  if numel(fields) > 1 && ~remove_field
	for i = 1:numel(fields)
	  if numel(fields{i}) == 2
		pth = M.(fields{i});
		text(pth(end,1), pth(end,2), fields{i},'color','w')
	  end
	end  
	pause(0.1);
	char = input('Do you want to add transects? y/n ', 's');
	if char == 'n'
	  fprintf('saving struct to %s.mat:\n', numel(fields), mask_name);
	  save([mask_name,'.mat'], '-struct', '-mat' , 'M');
	  return
	end
  end
  
  if ~remove_field
	figure(1);
	counter = 0;
	while true
	  counter = counter+1;

	  fprintf('\n  Set transect...\n');
	  [coords, status] = getcoordinates(tmp_maskp);
	  
	  if status == 1
		break;
	  end
	  
	  trpath = getpath(coords,int_maskp);

	  %% Display
	  mx= max(max(abs(maskp)));
	  for i = 1:size(trpath,1);
		tmp_maskp(trpath(i,2),trpath(i,1)) = mx+4*trpath(i,3);
	  end	
	  imagesc(tmp_maskp); set(gca,'ydir','normal');
	  drawnow
	  
	  disp(trpath)
	  fprintf('Transect %d\n',counter);
	  TName = input('Enter transect identifier: ','s');
	  
	  M.(TName) = trpath;
	end
  end
  
  if remove_field
	RField = input('Remove field: ', 's');
	M = rmfield(M,RField);
  end
  
  fields = fieldnames(M);
  fprintf('saving %d fields to %s.mat:\n', numel(fields), mask_name);
  for i = 1:numel(fields)
	fprintf('     %s\n', fields{i});
  end
  
  save([mask_name,'.mat'], '-struct', '-mat' , 'M');

end
