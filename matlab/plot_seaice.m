function [state,pars,add] = plot_seaice(fname, opts)

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

    M0   = 0;
    Mvar = 1.0;
    
    t0o  = 15;
    
    c1 = 3.8e-3;
    c2 = 21.87;
    c3 = 265.5;
    c4 = 17.67;
    c5 = 243.5;

    T0   = c3*c4*t0o / (c2*c5+(c2-c4)*t0o);
    Tvar = 1;

    scales = [Hvar, Qvar, Mvar, Tvar];
    backgr = [H0, Q0, M0, T0];
    
    [state,pars,add] = readhdf5(fname, si_nun, n, m, si_l, opts);

    titles = {'H','Q_T^{sa}','M','T'};

    simask = squeeze(state(3, :, :, :));
    
    for i = 1:si_nun
        figure(14+i)
        field = backgr(i) + scales(i)*squeeze(state(i, :, :, :));
        field(logical(surfm)) = NaN;
        
        diff = max(max(field))-min(min(field));
        fprintf('max(%s)-min(%s) = %f\n', titles{i}(1), titles{i}(1), ...
                diff);
        imagesc(RtD*x, RtD*(y), field'); hold on
        plot_mask(surfm,x, y); hold off


        set(gca, 'ydir', 'normal')
        title(titles{i});
        colorbar
    end
    
    if readFluxes
        plot_fluxes(add, 14+i+1);
    end
    
end