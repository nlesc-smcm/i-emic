function [state,pars,add,fluxes] = plot_seaice(fname, opts)

    if nargin < 2
        opts = [];
    end
    
    if nargin < 1
        fname = 'seaice_output.h5';
    end
    
    if isfield(opts, 'readFluxes')
        readFluxes = opts.readFluxes;
    else
        readFluxes = false;
    end

    if isfield(opts, 'readEV')
        readEV = opts.readEV;
    else
        readEV = false;
    end

    if isfield(opts, 'fig_ctr')
        fig_ctr = opts.fig_ctr;
    else
        fig_ctr = 21; % first figure handle number
    end

    [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');
    
    surfm      = landm(2:n+1,2:m+1,l+1);    % Only interior surface points
    landm_int  = landm(2:n+1,2:m+1,2:l+1);
    summask    = sum(landm_int,3);

    RtD = 180/pi;

    si_nun = 4;
    si_l   = 1;

    H0   = 0.01;
    Hvar = 1;
    Q0   = -100;
    Qvar = 498.8928;
    %Q0   = 0.0;
    % Qvar = 1.0;

    M0   = 0.0;
    Mvar = 1.0;
    
    t0o  = 15;
    t0i  = -15;
    
    c1 = 3.8e-3;
    c2 = 21.87;
    c3 = 265.5;
    c4 = 17.67;
    c5 = 243.5;

    T0   = t0i;
    Tvar = 1;

    scales = [Hvar, Qvar, Mvar, Tvar];
    backgr = [H0, Q0, M0, T0];
    
    [state,pars,add] = readhdf5(fname, si_nun, n, m, si_l, opts);

    titles = {'H','Q_T^{sa}','M','T'};

    simask = squeeze(state(3, :, :, :)) > 0.1;
    
    for i = 1:si_nun
        figure(fig_ctr); fig_ctr = fig_ctr+1;
        field = backgr(i) + scales(i)*squeeze(state(i, :, :, :));
        field(logical(surfm)) = NaN;

        if ~readEV
            field(~logical(simask)) = NaN;
        end
        
        diff = max(max(field))-min(min(field));
        fprintf('max(%s)-min(%s) = %f\n', titles{i}(1), titles{i}(1), ...
                diff);
        imagesc(RtD*x, RtD*(y), field'); hold on
        ca = caxis;
        plot_mask(surfm,x, y); hold off
        caxis(ca);


        set(gca, 'ydir', 'normal')
        title(titles{i});
        colorbar
    end

    opts.mask = logical(surfm) | ~logical(simask);
    fluxes = [];
    if readFluxes
        fluxes = plot_fluxes(add, fig_ctr, 'SeaIce: ', opts);
    end
    
end