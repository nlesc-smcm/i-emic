function [titles, cdata] = plot_cdata(fname, lsty)
    
    if nargin < 2
        lsty = 'k.-';
    end
    if (nargin < 1)
        fname = 'cdata.txt'
    end
    
    [titles, cdata] = load_cdata(fname);

    maxPsi = 0;
    minPsi = 0;
    
    for i = 2:size(cdata,2)
        figure(i)
        if strcmp(titles{i}, 'max(Psi)')
            maxPsi = i;
        elseif strcmp(titles{i}, 'min(Psi)')
            minPsi = i;
        end
        
        if strcmp(titles{i}, '||F||')
            semilogy(cdata(:,1),cdata(:,i),lsty);
            xlabel('par');
        else
            plot(cdata(:,1),cdata(:,i),lsty);
            xlabel('par');
        end
        title(titles{i})
        grid on;
    end
    
    if ((maxPsi > 0) && (minPsi > 0))
        figure(i+1);
        plot(cdata(:,1),cdata(:,maxPsi)+cdata(:,minPsi),lsty);
        title('max(Psi)+min(Psi)');
    end       

end