function [M] = read_many_forts(range, dim)
  obs = numel(range); % number of observations
  % data matrix
  M   = zeros(obs, dim);
  for k = 1:obs
	fprintf('reading fort.%d\n', range(k));
	[lab icp par xl xlp det sig sol solup soleig] = ...
    readfort3(0,strcat('fort.',num2str(range(k))), true);
	M(k,:) = sol(:);	
  end
  %M(:,4:6:end) = 0*M(:,4:6:end)/10000;
end
