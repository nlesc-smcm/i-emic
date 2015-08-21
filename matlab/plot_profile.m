function [profile] = plot_profile(filename, range)
  % we assume that you are in a directory below the main
  % project directory, for instance proj_dir/rundir
  system(['../scripts/gatherprofile.sh ' filename ' converted']);
  profile = importdata('converted');
  N       = profile.data(1); % number of experiments
  M       = numel(profile.textdata);
  cores   = 2.^(0:N-1);
  x_axis  = 1:N;
  if nargin == 1
	range = 1:M;
  end
  for i = range
	semilogy(x_axis, profile.data((i-1)*N+2:i*N+1),'.-','color',[rand rand rand]);
	hold on
  end
  hold off
  legend(profile.textdata(range));
	  
end
