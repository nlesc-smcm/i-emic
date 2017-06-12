function [sol, pars] = readhdf5(file, nun, n, m, l)

% load state
    sol = h5read(file, '/State/Values');
    sol = reshape(sol, nun, n, m, l);

    % load parameters
    info  = h5info(file, '/Parameters');
    npars = size(info.Datasets, 1);

    for i = 1:npars
        parname   = info.Datasets(i).Name;
        fieldname = regexprep(parname, ' |-', '_');
        pars.(fieldname) = ...
            h5read(file, ['/Parameters/' parname]);
    end
end