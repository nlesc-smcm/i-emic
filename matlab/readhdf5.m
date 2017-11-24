function [sol, pars, additional] = readhdf5(file, nun, n, m, l, opts)

    if nargin < 6
        opts.tmp = 0
    end
    
    if isfield(opts, 'readE')
        readE = opts.readE;
    else
        readE = false;
    end

    if isfield(opts, 'readP')
        readP = opts.readP;
    else
        readP = false;
    end
    
    if isfield(opts, 'readParameters')
        readPars = opts.readParameters;
    else
        readPars = false;
    end
    
    % read state
    sol = h5read(file, '/State/Values');
    
    dim = n*m*l*nun;
    
    sol = reshape(sol(1:dim), nun, n, m, l);
    
    % read parameters
    pars = [];
    if readPars
        
        info  = h5info(file, '/Parameters');
        npars = size(info.Datasets, 1);

        for i = 1:npars
            parname   = info.Datasets(i).Name;
            fieldname = regexprep(parname, ' |-', '_');
            pars.(fieldname) = ...
                h5read(file, ['/Parameters/' parname]);
        end
       
    end
    
    % read additional fields if requested
    additional = [];
    if readE
        additional.E = h5read(file, '/E/Values');
    end
    
    if readP
        additional.P = h5read(file, '/P/Values');
    end

end