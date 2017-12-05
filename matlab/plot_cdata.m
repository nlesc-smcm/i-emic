function [titles, cdata] = plot_cdata(fname, lsty)
    
    if nargin < 2
        lsty = 'k.-';
    end
    if (nargin < 1)
        fname = 'cdata.txt'
    end

    system(['tail -n +2 ', fname, ' > tmp']); 
    cdata = load('tmp');
    system(['head -n 1 ', fname, ' > tmp']);
    fid = fopen('tmp','r');

    ncol = size(cdata,2);

    titles = cell(ncol,1);

    for i = 1:ncol
        titles{i} = fscanf(fid,'%s',1);
    end

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