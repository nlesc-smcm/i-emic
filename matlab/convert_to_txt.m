function [] = convert_to_txt(mask_name, output_name)
  % convert Michiel's binary mask to a plaintext format
  % this scripts assumes that we never have more than 99 vertical
  % layers!

  M = load([mask_name, '.mat']);

  maskp = M.maskp;

  fid      =  fopen([output_name, '.txt'],'w');

  [m,n,l]  =  size(maskp);

  fprintf('  writing to %s\n',[output_name, '.txt']);
  fprintf(fid,'%3d %3d %3d\n',n,m,max(max(maskp)));
  for j = m:-1:1
	for i = 1:n
	  fprintf(fid,'%3d',maskp(j,i));
	end
	fprintf(fid, '\n');
  end
end
