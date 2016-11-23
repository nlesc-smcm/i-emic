function [] = add_basin(mask_name);
  %% add_basin(mask_name);
  %%
  %% Append a basin mask to mask file: mask_name
  %%  mask_name should be the name of a .mat file
  %%  with at least the field maskp
  %%
  %% Transect data is denoted with 2 capitals
  %% Basin data is identified with 4 capitals
		 
  M      = load([mask_name, '.mat']);
  maskp  = M.maskp;
  [m,n]  = size(maskp);
  tmp_maskp = maskp;
  int_maskp = maskp;
  
  fields = fieldnames(M);
  
  %% Display the basins already in there
  mx = max(max(abs(maskp)));
  for i = 1:numel(fields)
	if numel(fields{i}) == 4
	  fprintf('Found %s\n', fields{i});
	end
  end

  figure(1)
  imagesc(tmp_maskp);
  set(gca,'ydir','normal');
  set(gca,'ytick',linspace(1,m,5))
  grid on
  
  char = input('Do you want to add a basin? y/n ', 's');
  if char == 'n'
	return
  end
  
  counter = 0;
  while true


	fprintf('\n  Set basin...\n');
	[coords, status] = getcoordinates(tmp_maskp);
	
	if status == 1
	  break;
	end
	
	basinmask = getbasin(coords, int_maskp);
	imagesc(basinmask);
	set(gca,'ydir','normal');
	set(gca,'ytick',linspace(1,m,5))
	grid on
	
	char = input('Satsified? y/n/q ', 's');
	if char == 'n'
	  continue
	elseif char == 'q'
	  break
	end	

	counter = counter+1;
	fprintf('Basin %d\n',counter);
	BName = input('Enter identifier (4 capitals NOAT,SOAT,...): ','s');
	
	M.(BName) = basinmask;
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
