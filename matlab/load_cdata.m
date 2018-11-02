function [titles, cdata] = load_cdata(fname)

    if (nargin < 1)
        fname = 'cdata.txt'
    end

    system(['tail -n +2 ', fname, ' > tmp']); 
    cdata = load('tmp');
    system(['head -n 1 ', fname, ' > tmp']);
    fid = fopen('tmp','r');

    ncol = size(cdata,2);

    titles = cell(ncol,1);

    i = 1;
    while i <= ncol
        titles{i} = fscanf(fid,'%s',1);

        % ignore initial pound
        if strcmp(titles{i},'#')
            i = i - 1
        end
        i = i + 1;                
    end
end