function [] = rename_transect(mask_name, old_name, new_name);
		 
  M = load([mask_name, '.mat']);

  fprintf('renaming %s to %s...\n', old_name, new_name);
  %% swap data
  data = M.(old_name);
  M.(new_name) = data;
  M = rmfield(M, old_name);

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
