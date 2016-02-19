function [profile] = plot_profile(filename, range, domain)
  % we assume that you are in a directory below the main
  % project directory, for instance proj_dir/rundir
  system(['../scripts/gatherprofile.sh ' filename ' converted']);
  profile = importdata('converted');
  N       = profile.data(1); % number of experiments
  M       = numel(profile.textdata);
  cores   = 2.^(1:N);
  x_axis  = 1:N;
%   x_axis  = cores;
  
  if nargin < 3
      domain = 1:N;
      if nargin < 2
          range = 1:M;
      end
  end 
  
  colm = lines(numel(range));
  ctr  = 1;
  linst = {'.-','.:','.--'};
  for i = range
      profile.textdata{i} = ['(',num2str(i),'): ',profile.textdata{i}];
      
      pllt = profile.data((i-1)*N+2:i*N+1);
      
      figure(8) 
      plot(x_axis(domain), pllt(domain),linst{mod(ctr,numel(linst))+1},...
                    'linewidth',2,'markersize',15,'color',colm(ctr,:));
      x = x_axis(domain(end));    
      y = pllt(domain(end));
      hold on
      text((1.02+(rand-.5)/30)*x,y,num2str(i),'color',colm(ctr,:));
      hold on
      
      figure(9) 
      semilogy(x_axis(domain), pllt(domain),linst{mod(ctr,numel(linst))+1},...
                    'linewidth',2,'markersize',15,'color',colm(ctr,:));
      x = x_axis(domain(end));    
      y = pllt(domain(end));
      hold on
      text((1.02+(rand-.5)/30)*x,y,profile.textdata{i},'color',colm(ctr,:));
      xlim([0,1.4*x_axis(end)]);
      hold on
      
      figure(10)
      plot(x_axis(domain),pllt(1)./pllt(domain),linst{mod(ctr,numel(linst))+1},...
                    'linewidth',2,'markersize',15,'color',colm(ctr,:));
      x = x_axis(domain(end));
      y = pllt(1)./pllt(domain(end));      
      hold on
      text((1.02+(rand-.5)/30)*x,y,num2str(i),'color',colm(ctr,:));
      
      ctr = ctr + 1;
  end
  figure(8)
  hold off
  figure(9)
  hold off
  figure(10) 
  legend(profile.textdata(range),'location','northwest');
  plot(x_axis(domain), x_axis(domain), 'k--')
  hold off
 
    
end
