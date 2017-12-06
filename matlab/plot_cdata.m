function [titles, cdata] = plot_cdata(fname, lsty)
    
    if nargin < 2
        lsty = 'k.-';
    end
    if (nargin < 1)
        fname = 'cdata.txt'
    end
    
    [titles, cdata] = load_cdata(fname);

    for i = 2:size(cdata,2)
        figure(i)
        if strcmp(titles{i}, 'NR') || strcmp(titles{i}, 'ds') || ...
                strcmp(titles{i}, '||F||') || strcmp(titles{i}, 'MV')
            plot(cdata(:,1),cdata(:,i),lsty);
            xlabel('par');
        else
            plot(cdata(:,1),cdata(:,i),lsty);
            xlabel('par');
        end
        title(titles{i})
        grid on;
    end

end