function [profile] = plot_profile(filename, range, domain)
  % we assume that you are in a directory below the main
  % project directory, for instance proj_dir/rundir
  system(['../scripts/gatherprofile.sh ' filename ' converted']);
  profile = importdata('converted');
  N       = profile.data(1) % number of experiments
  M       = numel(profile.textdata)
  cores   = 2.^(0:N-1);
  x_axis  = 1:N;
  x_axis  = cores;
  
  if nargin < 3
      domain = 1:N;
      if nargin < 2
          range = 1:M;
      end
  end 
  
  colm = lines(numel(range));
  ctr  = 1;
  for i = range
      pllt = profile.data((i-1)*N+2:i*N+1);
      figure(9) 
      semilogy(x_axis(domain), pllt(domain),'.-','linewidth',2,'markersize',15,'color',colm(ctr,:));
      hold on
      
      figure(10)
      plot(x_axis(domain),pllt(1)./pllt(domain),'.-','linewidth',2,'markersize',15,'color',colm(ctr,:));
      hold on
      
      ctr = ctr + 1;
  end
  figure(9)
  hold off
  figure(10) 
  legend(profile.textdata(range),'location','northwest');
  plot(x_axis(domain), x_axis(domain), 'k--')
  hold off
 
    
end
