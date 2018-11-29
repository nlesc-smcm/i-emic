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
    
    if isfield(opts, 'sifield')
        sifield  = opts.sifield;
        plot_all = false;
    else
        sifield = '';
        plot_all = true;
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
    
    if isfield(opts,'exporteps')
        exporteps = opts.exporteps;
    else
        exporteps = false;
    end

    if isfield(opts,'invert')
        invert = opts.invert;
    else
        invert = false;
    end
    
    if isfield(opts,'input_caxis')
        set_caxis = true;
    else
        set_caxis = false;
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

    simask = squeeze(state(3, :, :, :)) > 0.2;
    
    for i = 1:si_nun
        
        if plot_all || strcmp(titles{i}, sifield)
            
            figure(fig_ctr); fig_ctr = fig_ctr+1;
            set(gca,'color',[.5,.5,.5]);
            
            field = backgr(i) + scales(i)*squeeze(state(i, :, :, :));
            field(logical(surfm)) = NaN;

            if ~readEV
                field(~logical(simask)) = NaN;
            end
            
            diff = max(max(field))-min(min(field));
            fprintf('max(%s)-min(%s) = %f\n', titles{i}(1), titles{i}(1), diff);
            %            imagesc(RtD*x, RtD*(y), field','AlphaData',~isnan(field')); ...
            %   hold on
            imagesc(RtD*x, RtD*(y), field'); ...
                hold on
            
            cnan = [100 100 100]/256;
            cbeg = [240 240 240]/256;
            cmid = [150 180 200]/256;
            cend = [180	220	240]/256;
            cmap  = [cnan; ...
                    linspace(cbeg(1),cmid(1),128)',...
                    linspace(cbeg(2),cmid(2),128)',...
                    linspace(cbeg(3),cmid(3),128)';...
                    linspace(cmid(1),cend(1),128)',...
                    linspace(cmid(2),cend(2),128)',...
                    linspace(cmid(3),cend(3),128)'];
            
            colormap(cmap);
            
            amin=min(caxis);
            amax=max(caxis);
            n = size(cmap,1);

            dmap=2*(amax-amin)/n;
            caxis([amin-dmap,amax]);
            
            %colormap(my_colmap(caxis,mean(caxis),128,[0,0.4470,0.7410],[0.7,0.7,0.7]))
            if set_caxis
                caxis(opts.input_caxis)
            end
            
            ca = caxis;
            plot_mask(surfm,x, y); 
            
            hold off
            caxis(ca);
            
            set(gca, 'ydir', 'normal')
            title(titles{i});
            xlabel('Longitude')
            ylabel('Latitude')
            colorbar
            
            if ~exporteps && invert
                invertcolors()
            end
            
            if exporteps
                exportfig(['seaice',titles{i},'.eps'],10,[19,11], ...
                          invert)
            end
        end
    end

    opts.mask = logical(surfm) | ~logical(simask);
    fluxes = [];
    if readFluxes
        fluxes = plot_fluxes(add, fig_ctr, 'SeaIce: ', opts);
    end
    
end