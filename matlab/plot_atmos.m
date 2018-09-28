function [state,pars,add] = plot_atmos(fname, opts)

    if nargin < 2
        opts = [];
        plot_default = true;
    else
        plot_default = false;
    end
    
    if nargin < 1
        fname = 'atmos_output.h5';
    end
    
    if nargin < 2 % defaults
        readEP  = true;
        readLST = true;      
        
        opts.readEP  = readEP;
        opts.readLST = readLST;
        opts.lst     = readLST;
    end
    
    if isfield(opts, 'EmP') 
        plotEmP = opts.EmP;
    else
        plotEmP = false;
    end
    opts.readEP = plotEmP;

    % land surface temperature
    if isfield(opts, 'lst') 
        readLST = opts.lst;
        opts.readLST = readLST;
    else
        readLST = false;
    end

    % sea ice surface temperature
    if isfield(opts, 'sit') 
        readSIT = opts.sit;
        opts.readSIT = readSIT;
    else
        readSIT = false;
    end

    if isfield(opts, 'readFluxes')
        readFluxes = opts.readFluxes;
    else
        readFluxes = false;
    end

    if isfield(opts, 'Ta')
        plot_ta = opts.Ta;
    else
        plot_ta = false;
    end

    if isfield(opts, 'humidity')
        plot_humidity = opts.humidity;
    else
        plot_humidity = false;
    end

    if isfield(opts, 'albedo')
        plot_albedo = opts.albedo;
    else
        plot_albedo = false;
    end

    if isfield(opts, 'invert')
        invert = opts.invert;
    else
        invert = false;
        opts.invert = invert;
    end
    
    if isfield(opts, 'fig_ctr')
        fig_ctr = opts.fig_ctr;
    else
        fig_ctr = 11; % first figure handle number
    end

    if isfield(opts, 'qsat')
        plot_qsat = opts.qsat;
        readLST   = opts.qsat;
        opts.readLST = readLST;
    else
        plot_qsat = false;
    end

    if isfield(opts, 'exportfig')
        export_to_file = opts.exportfig;
    else
        export_to_file = false;
    end

    [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');
    surfm      = landm(2:n+1,2:m+1,l+1);    % Only interior surface points
    landm_int  = landm(2:n+1,2:m+1,2:l+1);
    summask = sum(landm_int,3);

    srf = [];
    greyness = .85;
    srf(:,:,1) = (1-greyness*(surfm'));
    srf(:,:,2) = (1-greyness*(surfm'));
    srf(:,:,3) = (1-greyness*(surfm'));

    atmos_nun = 3;
    atmos_l   = 1;
    
    [state, pars, add] = readhdf5(fname, atmos_nun, n, m, atmos_l,opts);

    if plotEmP
        E = reshape(add.E,n,m);
        P = reshape(add.P,n,m);
    end

    if readLST
        LST = reshape(add.LST, n, m);
        SST = reshape(add.SST, n, m);
    end

    if readSIT
        SIT = reshape(add.SIT, n, m);
        MSI = reshape(add.MSI, n, m);
    else
        SIT = 0;
        MSI = 0;
    end

    RtD = 180/pi;
    
    % reference values
    T0o  = 15.0;
    T0i  = -5;
    q0   = 8e-3;
    a0   = 0.0;
    
    % typical deviations
    tdim = 1;
    qdim = 1e-3;
    adim = 1;

    % constants
    rhoa = 1.25;
    rhoo = 1024;
    ce   = 1.3e-03;
    uw   = 8.5;
    eta  = (rhoa / rhoo) * ce * uw;
    year = 3600*24*365;
    
    fprintf('  q0 = %e (kg / kg)\n', q0);        
    fprintf('qdim = %e (kg / kg)\n', qdim);
    fprintf(' eta = %e \n', eta);

    c1 = 3.8e-3;
    c2 = 21.87;
    c3 = 265.5;
    c4 = 17.67;
    c5 = 243.5;
    
    Ta  = T0o + tdim * squeeze(state(1,:,:,:));
    qa  = q0 + qdim * squeeze(state(2,:,:,:));
    %qa  =  qdim * squeeze(state(2,:,:,:));
    qa(logical(surfm)) = NaN; % humidity above land points does not contribute 
    Aa  = a0 + adim * squeeze(state(3,:,:,:));
    Tz  = mean(Ta,1); % zonal mean
    qz  = mean(qa,1); % zonal mean

    if plot_ta || plot_default
        figure(fig_ctr); fig_ctr = fig_ctr+1;

        img = Ta';
        Tdiff = max(max(Ta))-min(min(Ta));
        fprintf('max(T) - min(T) = %f\n', Tdiff);
        
        %imagesc(RtD*x,RtD*(y),img); hold on
        %contourf(RtD*x,RtD*(y),img,10,'linestyle','none','Fill','on'); hold on;
        %plot_mask(surfm,x,y); hold on
        %contour(RtD*x,RtD*(y),img,20,'linecolor','k',...
        %        'linewidth',1.2,'linestyle','-'); hold on;
        
        contourf(RtD*x,RtD*(y),img,15,... 
                 'linewidth',1.0,'linestyle','-'); hold on;
        contour(RtD*x,RtD*(y),surfm',1,'linecolor','k','linewidth',2.0,'linestyle','-'); hold on;
                
        set(gca,'ydir','normal')
        colorbar
        cmap = my_colmap(caxis);
        colormap(cmap)
        colorbar

        hold off
        
        title('Atmospheric temperature')
        xlabel('Longitude')
        ylabel('Latitude')
        
        if export_to_file
            exportfig('atmosTemp.eps',10,[14,10],invert)
        end
    end

    if plot_humidity || plot_default

        figure(fig_ctr); fig_ctr = fig_ctr+1;
        img = qa';

        imagesc(RtD*x,RtD*(y),img); hold on                
        c = contour(RtD*x,RtD*(y),img,15,'Visible', 'on', ...
                    'linewidth',1.0,'linecolor','k'); hold on
        plot_mask(surfm,x,y); hold on

        hold off
        set(gca,'ydir','normal')
        
        colorbar
        cmap = my_colmap(caxis);
        colormap(cmap)
        
        drawnow
        title('Humidity anomaly (kg / kg)')
        xlabel('Longitude')
        ylabel('Latitude')
        
        if export_to_file
            exportfig('atmosq.eps',10,[14,10],invert)
        end
    end
    
    if plot_albedo || plot_default
        figure(fig_ctr); fig_ctr = fig_ctr+1;
        img = Aa';
        imagesc(RtD*x,RtD*(y),img); 
        set(gca,'ydir','normal')
        
        %cmap = my_colmap(caxis);
        colormap(gray)
        colorbar

        title('Albedo')
        xlabel('Longitude')
        ylabel('Latitude')
        
        if export_to_file
            exportfig('atmosA.eps',10,[14,10],invert)
        end
    end
    
    if plotEmP

        figure(fig_ctr); fig_ctr = fig_ctr+1;        

        Pmy = P*3600*24*365;
        img = Pmy';
        img(img == 0) = NaN;

        %contour(RtD*x,RtD*(y),surfm',1,'linecolor','k','linewidth',2.0,'linestyle','-'); hold on;
        c = contourf(RtD*x,RtD*(y),img,12,'k','Visible', 'on', ...
                     'linewidth',.5); 
                
        set(gca,'ydir','normal')
        cmap = my_colmap(caxis);
        colormap(cmap)
        colorbar
        
        title('P (m/y)')
        xlabel('Longitude')
        ylabel('Latitude')

        figure(fig_ctr); fig_ctr = fig_ctr+1;        

        Emy = E*3600*24*365;
        img = Emy';
        img(img == 0) = NaN;
                
        %contour(RtD*x,RtD*(y),surfm',1,'linecolor','k','linewidth',2.0,'linestyle','-'); hold on;
        c = contourf(RtD*x,RtD*(y),img,12,'k','Visible', 'on', ...
                     'linewidth',.5);
                
        set(gca,'ydir','normal')
        cmap = my_colmap(caxis);
        colormap(cmap)
        colorbar
        
        title('E (m/y)')
        xlabel('Longitude')
        ylabel('Latitude')

        if export_to_file
            exportfig('atmosE.eps',10,[14,10],invert)
        end

        EmPmy = (E-P)*3600*24*365;
        img   = EmPmy';
        img(img == 0) = NaN;
        
        figure(fig_ctr); fig_ctr = fig_ctr+1;        
        %imagesc(RtD*x,RtD*(y),img); 

         c = contourf(RtD*x,RtD*(y),img,12,'k','Visible', 'on', ...
                            'linewidth',.5);
        % hold off
        set(gca,'ydir','normal')
        cmap = my_colmap(caxis,0);
        colormap(cmap)
        colorbar

        title('E-P (m/y)')
        xlabel('Longitude')
        ylabel('Latitude')
        
        if export_to_file
            exportfig('atmosEmP.eps',10,[14,10],invert)
        end
    end
    
    if readLST || plot_default
        figure(fig_ctr); fig_ctr = fig_ctr+1;
        imagesc(RtD*x,RtD*(y), LST' + T0o);
        set(gca,'ydir', 'normal');
        cmap = my_colmap(caxis);
        colormap(cmap)
        colorbar
        title('Surface temperature')
        xlabel('Longitude')
        ylabel('Latitude')
    end
    
    if plot_qsat
        figure(fig_ctr); fig_ctr = fig_ctr+1;
        sst  = SST' + T0o;
        qso  = @(T) c1 * exp(c4 * T ./ (T + c5));
        qsi  = @(T) c1 * exp(c2 * T ./ (T + c3));
        
        dqso     = (qso(T0o+1e-6)-qso(T0o))/1e-6;
        dqso     = 5e-4;
        dqsi     = (qsi(T0i+1e-6)-qsi(T0i))/1e-6;
        qso_lin  = qso(T0o) + dqso * SST';
        qsi_lin  = qsi(T0i) + dqsi * SIT';

        fprintf(' background qsat over ocean   = %2.1f g/kg\n', qso(T0o)*1e3);
        fprintf(' background qsat over sea ice = %2.1f g/kg\n\n', qsi(T0i)*1e3);
        fprintf(' d/dT qsat over ocean   = %1.3e\n', dqso);
        fprintf(' d/dT qsat over sea ice = %1.3e\n\n', dqsi);
        
        qs = qso_lin .* ( 1-MSI') + MSI' .* qsi_lin;
       
        % full nonlinear qso
        %qso_nonlin  = qso(sst);
        
        contourf(RtD*x,RtD*(y),qs,12); hold off
        title('q_{sat} (kg/kg)', 'interpreter', 'TeX');
        xlabel('Longitude');
        ylabel('Latitude');
        cmap = [my_colmap(caxis)];
        colormap(cmap)
        colorbar
    end
    
    if plot_qsat && plot_humidity
        figure(fig_ctr); fig_ctr = fig_ctr+1;
        contourf(RtD*x,RtD*(y),eta*(qs-qa')*year,12); hold off
        title('q_{sat} - q', 'interpreter', 'TeX');
        xlabel('Longitude');
        ylabel('Latitude');
        cmap = [my_colmap(caxis)];
        colormap(cmap)
        colorbar

        fprintf(' background evaporation Eo0 = %e m/s = %2.1f m/y\n', ...
                eta * ( qso(T0o) - q0), ...
                eta * ( qso(T0o) - q0)*year);
        fprintf(' background sublimation Eoi = %e m/s = %2.1f m/y\n', ...
                eta * ( qsi(T0i) - q0),...
                eta * ( qsi(T0i) - q0)*year);
    end
    
    if readFluxes
        out = plot_fluxes(add, fig_ctr,'Atmos: ', opts); 
    end
    
end