function [] = convert_to_bin(mask_name)

% convert plaintext flat mask to a binary format
% this scripts assumes that we never have more than 99 vertical
% layers!

fid  =  fopen([mask_name, '.txt'],'r');
fprintf('  opening %s\n',[mask_name, '.txt']);

n = fscanf(fid,'%3d',1)
m = fscanf(fid,'%3d',1)
levels = fscanf(fid,'%3d',1)

maskp = zeros(m,n);

for j = m:-1:1
  for i = 1:n
	maskp(j,i) = fscanf(fid,'%3d',1);
  end
end

imagesc(maskp)
set(gca,'ydir','normal');


fprintf('  writing to %s\n',[mask_name, '.mat']);
save([mask_name,'.mat'],'maskp');
