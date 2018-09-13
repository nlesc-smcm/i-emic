function [sol, pars, additional] = readhdf5(file, nun, n, m, l, opts)
    
    if nargin < 6
        opts.tmp = 0;
    end
    
    if nargin < 2
        no_reshape = true;
    else
        no_reshape = false;
    end
    
    if isfield(opts, 'readEP')
        readEP = opts.readEP;
    else
        readEP = false;
    end
    
    if isfield(opts, 'readEV')
        readEV = opts.readEV;
    else
        readEV = false;
    end

    if isfield(opts, 'readLST')
        readLST = opts.readLST;
    else
        readLST = false;
    end

    if isfield(opts, 'readFluxes')
        readFluxes = opts.readFluxes;
    else
        readFluxes = false;
    end

    if isfield(opts, 'readParameters')
        readPars = opts.readParameters;
    else
        readPars = false;
    end

    % initialize additional fields if requested
    additional = [];
    
    %------------------------------------------------------------------
    %------------------------------------------------------------------    
    % read state
    
    if readEV % read eigenvector and values
        kmax = h5read(file, '/MetaData/NumEigs');

        if isfield(opts, 'evindex')
            evindex = opts.evindex;            
        else
            evindex = 0;
        end
        
        if (evindex >= kmax || evindex < 0)
            fprintf(['WARNING: evindex not available, reading ev 0.\' ...
                     'n']);
            evindex = 0;
        end

        sol = h5read(file, ['/EV_Real_',num2str(evindex),['/' ...
                            'Values']]);

        alphaRe = h5read(file, '/EigenValues/AlphaRe');
        alphaIm = h5read(file, '/EigenValues/AlphaIm');
        betaRe  = h5read(file, '/EigenValues/BetaRe' );
        betaIm  = h5read(file, '/EigenValues/BetaIm' );

        additional.alpha = alphaRe + alphaIm*1i;
        additional.beta  = betaRe + betaIm*1i;
        
    else % read normal state
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
    
    if readEP
        additional.E = h5read(file, '/E/Values');
        additional.P = h5read(file, '/P/Values');
    end
    
    if readLST
        additional.LST = h5read(file, '/lst/Values');
        additional.SST = h5read(file, '/sst/Values');
    end

    if readFluxes
        info = h5info(file);

        % check for groups containing the characters 'Flux'
        % and add the contents to the struct
        nGroups = numel(info.Groups);
        for i = 1:nGroups
            groupName = info.Groups(i).Name;
            match = regexp(groupName, '.*Flux', 'match');
            if numel(match) > 0
                field  = match{1}(2:end);
                values = h5read(file, ['/',field,'/Values']);
                additional = setfield(additional, field, values);
            end
        end
    end
end