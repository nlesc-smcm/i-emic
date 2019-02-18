function [profile] = plot_profile(filename, range, start)

    system(['/home/erik/Projects/i-emic/scripts/gatherprofile.sh ' filename ' converted']);
    profile = importdata('converted');

    N       = profile.data(1); % number of experiments
    M       = numel(profile.textdata);
    cores   = [4,8,16,24];
    x_axis  = 1:N;
    x_axis  = cores;

    if nargin < 3
        domain = 1:N;
        start = 1;
        if nargin < 2
            range = 1:M;
        end
    else
        domain = start:N;
    end
    
    
    colm = lines(numel(range));
    ctr  = 1;
    linst = {'.-','.:','.--'};
    for i = range
        profile.textdata{i} = ['(',num2str(i),'): ',profile.textdata{i}];
        
        pllt = profile.data((i-1)*N+2:i*N+1);
        
        %       figure(8) 
        %       plot(x_axis(domain), pllt(domain),linst{mod(ctr,numel(linst))+1},...
        %                     'linewidth',2,'markersize',15,'color',colm(ctr,:));
        %       x = x_axis(domain(end));    
        %       y = pllt(domain(end));
        %       hold on
        %       text((1.02+(rand-.5)/30)*x,y,num2str(i),'color',colm(ctr,:));
        %       hold on
        
        figure(9) 
        semilogy(x_axis(domain), pllt(domain),linst{mod(ctr,numel(linst))+1},...
               'linewidth',2,'markersize',15,'color',colm(ctr,:));
        x = x_axis(domain(end));    
        y = pllt(domain(end));
        hold on
        text((1.02+(rand-.5)/30)*x,y,[profile.textdata{i},': ',num2str(y)],'color',colm(ctr,:));
        xlim([0,1.4*x_axis(end)]);
        hold on
        
        figure(10)
        plot(x_axis(domain),pllt(start)./pllt(domain), linst{mod(ctr,numel(linst))+1},...
             'linewidth',2,'markersize',15,'color',colm(ctr,:));
        x = x_axis(domain(end));
        y = pllt(1)./pllt(domain(end));      
        hold on
        text((1.02+(rand-.5)/30)*x,y,num2str(i),'color',colm(ctr,:));
        
        ctr = ctr + 1;
    end
    %   figure(8)
    %   hold off
    figure(9)
    hold off
    xlabel('#procs');
    ylabel('time');
    exportfig('profile1.eps',10,[24,12]);
    
    figure(10) 

    plot(x_axis(domain), x_axis(domain)./x_axis(start), 'k--')
    legend(profile.textdata(range),'location','northwest');
    
    hold off
    xlabel('#procs');
    ylabel('speedup');
    exportfig('profile2.eps',10,[24,12]);
    
end
