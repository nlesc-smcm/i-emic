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

    for i = 1:ncol
        titles{i} = fscanf(fid,'%s',1);
    end    
end