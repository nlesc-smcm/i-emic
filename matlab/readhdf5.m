function [sol, pars, additional] = readhdf5(file, nun, n, m, l, opts)
    
    if nargin < 6
        opts.tmp = 0
    end
    
    if nargin < 2
        no_reshape = true;
    else
        no_reshape = false;
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
    
    if isfield(opts, 'readEV')
        readEV = opts.readEV;
    else
        readEV = false;
    end

    if isfield(opts, 'salflux')
        readSalFlux = opts.salflux;
    else
        readSalFlux = false;
    end

    if isfield(opts, 'temflux')
        readTemFlux = opts.temflux;
    else
        readTemFlux = false;
    end

    if isfield(opts, 'readParameters')
        readPars = opts.readParameters;
    else
        readPars = false;
    end
    
    %------------------------------------------------------------------
    %------------------------------------------------------------------    
    % read state
    
    if readEV
        kmax = 4 % hardcoded for now

        if isfield(opts, 'eignum')
            eign = opts.eignum;
        else
            eign = 0;
        end

        sol = h5read(file, ['/EV_Real_',num2str(eign),'/Values']);
    else
        sol = h5read(file, '/State/Values');
    end
    
    if ~no_reshape
        dim = n*m*l*nun;
        sol = reshape(sol(1:dim), nun, n, m, l);
    end
    
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
    
    if readSalFlux
        additional.SalFlux = h5read(file, '/SalinityFlux/Values');
    end

    if readTemFlux
        additional.TemFlux = h5read(file, '/TemperatureFlux/Values');
    end    
end