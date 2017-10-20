function [state,pars,add] = plot_atmos(fname, opts)

    if nargin < 2
        opts = [];
    end

    if nargin < 1
        fname = 'atmos_output.h5';

        readE = true;
        readP = true;      
        opts.readE = readE;
        opts.readP = readP;
    end
    
    if isfield(opts, 'readE')
        readE = opts.readE;
    else
        readE = false;
        opts.readE = readE;
    end

    if isfield(opts, 'readP')
        readP = opts.readP;
    else
        readP = false;
        opts.readP = readP;
    end

    [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');
    surfm      = landm(2:n+1,2:m+1,l+1);    % Only interior surface points
    landm_int  = landm(2:n+1,2:m+1,2:l+1);

    srf = [];
    greyness = .85;
    srf(:,:,1) = (1-greyness*(surfm'));
    srf(:,:,2) = (1-greyness*(surfm'));
    srf(:,:,3) = (1-greyness*(surfm'));

    atmos_nun = 2;
    atmos_l = 1;
    
    [state,pars,add] = readhdf5(fname, atmos_nun, n, m, atmos_l,opts);
    
    if readE
        E = reshape(add.E,n,m);
    end
    
    if readP
        P = reshape(add.P,n,m);
    end
    
    T0  = 15.0;   %//! reference temperature
    q0  = 0.0015;
    RtD = 180/pi;    
    
    tdim = 1;
    qdim = 0.01;
    
    Ta  = T0 + tdim * squeeze(state(1,:,:,:));
    qa  = q0 + qdim * squeeze(state(2,:,:,:));
    Tz  = mean(Ta,1); % zonal mean
    qz  = mean(qa,1); % zonal mean

    figure(9)

    img = Ta';
    contourf(RtD*x,RtD*(y),img,20,'Visible','off'); hold on;
    image(RtD*x,RtD*(y),srf,'AlphaData',.2);
    c = contour(RtD*x,RtD*(y),img,20,'Visible', 'on','linewidth',1);
    colorbar
    caxis([min(min(Ta)),max(max(Ta))])
    hold off
    drawnow
    title('Atmospheric temperature')
    xlabel('Longitude')
    ylabel('Latitude')
    exportfig('atmosTemp.eps')

    figure(10)
    img = (Ta-repmat(Tz,n,1))';
    contourf(RtD*x,RtD*(y),img,20,'Visible','off'); hold on;
    image(RtD*x,RtD*(y),srf,'AlphaData',.2);
    c = contour(RtD*x,RtD*(y),img,20,'Visible', 'on','linewidth',1);
    colorbar
    caxis([min(min(img)),max(max(img))])
    hold off
    drawnow
    title('Ta anomaly')
    xlabel('Longitude')
    ylabel('Latitude')

    figure(11)
    img = (qa-repmat(qz,n,1))';
    contourf(RtD*x,RtD*(y),img,20,'Visible','off'); hold on;
    image(RtD*x,RtD*(y),srf,'AlphaData',.2);
    c = contour(RtD*x,RtD*(y),img,20,'Visible', 'on','linewidth',1);
    colorbar
    caxis([min(min(img)),max(max(img))])
    hold off
    drawnow
    title('qa anomaly')
    xlabel('Longitude')
    ylabel('Latitude')
    
    figure(12)
    img = qa';
    contourf(RtD*x,RtD*(y),img,20,'Visible','off'); hold on;
    image(RtD*x,RtD*(y),srf,'AlphaData',.2);
    c = contour(RtD*x,RtD*(y),img,20,'Visible', 'on','linewidth',1);
    colorbar
    caxis([min(min(qa)),max(max(qa))])
    hold off
    drawnow
    title('Atmospheric humidity (kg / kg)')
    xlabel('Longitude')
    ylabel('Latitude')
    exportfig('atmosq.eps')
    
    if readE && readP
        figure(13) 
        EmP = E-P;
        EmP(EmP==0)=NaN;
        img = EmP';
        contourf(RtD*x,RtD*(y),img,10,'Visible','off'); hold on;
        image(RtD*x,RtD*(y),srf,'AlphaData',.2);
        c = contour(RtD*x,RtD*(y),img,20,'Visible', 'on','linewidth',1);
        colorbar

        hold off
        drawnow
        title('E-P')
        xlabel('Longitude')
        ylabel('Latitude')
    end
    
end