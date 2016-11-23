function [] = remove_field(mask_name, field_name, force)

  if nargin < 3
	force = false;
  end

  fprintf('Loading %s\n', mask_name);
  M = load([mask_name, '.mat']);
  
  prompt = sprintf('Do you really want to remove %s? y/n\n',field_name);
  
  char = input(prompt, 's');
  if char == 'y'
	fprintf('Removing %s from %s\n', field_name, mask_name);
	M = rmfield(M, field_name);
	fprintf(['Saving to ', mask_name,'.mat\n']);
	onOctave = (exist ('OCTAVE_VERSION', 'builtin') > 0);
	
	if onOctave	  
	  save([mask_name,'.mat'], '-struct', '-mat' , 'M');
	else
	  save([mask_name,'.mat'], '-struct', 'M');
	end
  end
  
end
