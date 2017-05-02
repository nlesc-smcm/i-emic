function [] = plot_ocean(solfile, maskfile, opts)
%---------------------------------------------------------------------
% PLOTTHCM - Mother script for plotting THCM output
%  usage: plot_ocean(solfile, maskfile, opts)
%
%  Father is M. den Toom, who conceived it 06-11-08
%  Messed up by Erik, 2015/2016/2017 -> t.e.mulder@uu.nl
%---------------------------------------------------------------------

    if nargin < 1
        solfile = 'fort.3';
    end

    if nargin < 2
        maskfile = 'fort.44';
        specify_mask = false;
    else
        specify_mask = true;
    end

    if nargin < 3
        opts.everything = true;
    end

    if isfield(opts, 'title_add')
        plot_title = true;
    else
        plot_title = false;
    end

    if isfield(opts, 'fname_add')
        export_to_file = true;
    else
        export_to_file = false;
    end

    % interpolation mode
    if isfield(opts, 'solfile2') && isfield(opts, 'maskfile2') ...
        && isfield(opts, 'interp_par')
        interp_mode = true;
        k           = opts.interp_par;
        solfile2    = opts.solfile2;
        maskfile2   = opts.maskfile2;
        fprintf(' --- interp mode enabled\n');
    else
        interp_mode = false;
    end

    if isfield(opts, 'fix_caxis')
        fix_caxis = true;
        fprintf(' --- fixing caxis\n');
    else
        fix_caxis = false;
    end

    fprintf(1,'----------------------------------------------\n')

    % - DEFINE CONSTANTS - ----------------------------------------------

    udim  = 0.1;                 %[m/s]    Velocity scale
    r0dim = 6.4e6;               %[m]      Radius of Earth
    T0    = 15;                  %[deg C]  Reference temperature
    S0    = 35;                  %[psu]    Reference salinity
    RtD   = 180/pi;              %[-]      Radians to degrees

    % - READ MASK - -----------------------------------------------------

    [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = ...
        readfort44(maskfile);

    if interp_mode
        [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm2] = ...
            readfort44(maskfile2);
        landm = (1-k)*landm + k*landm2;
    end

    surfm      = landm(2:n+1,2:m+1,l+1);  %Only interior surface points
    landm_int  = landm(2:n+1,2:m+1,2:l+1);
    dx         = (xu(n+1)-xu(1))/n;
    dy         = (yv(m+1)-yv(1))/m;
    dz         = (zw(l+1)-zw(1))/l;

    % - Create surface landmask image
    summask = sum(landm_int,3);
    summask = summask / max(max(abs(summask)));
    summask = summask.^3;

    % - Deduce grid stretching
    [qz,dfzt,dfzw] = gridstretch(zw);

    % - READ SOLUTION
    if strcmp(solfile(end-1:end),'h5')       % (.h5 version)
        [sol] = readhdf5(solfile, nun, n, m, l);
    else % (fort.3 version)
        [~,~,~,~,~,~,~,sol,~,~] = readfort3(la, solfile);
    end

    if interp_mode
        if strcmp(solfile2(end-1:end),'h5')       % (.h5 version)
            [sol2] = readhdf5(solfile2, nun, n, m, l);
        else % (fort.3 version)
            [~,~,~,~,~,~,~,sol2,~,~] = readfort3(la, solfile2);
        end
        sol = (1-k)*sol+k*sol2;
    end

    % - EXTRACT SOLUTION COMPONENTS - -----------------------------------
    [u,v,w,p,T,S] = extractsol(sol);

    % --- Create colormaps

    par = [0    0.4470    0.7410;  0.8500    0.3250    0.0980];
    neg  = par(1,:);
    pos  = par(2,:);
    N1   = 64;
    N2   = 64;
    mid  = [1,1,1];
    col1 = [linspace(neg(1),mid(1),N1)',linspace(neg(2),mid(2),N1)',linspace(neg(3),mid(3),N1)'];
    col2 = [linspace(mid(1),pos(1),N2)',linspace(mid(2),pos(2),N2)',linspace(mid(3),pos(3),N2)'];
    col_white = [col1;col2];

    mid  = [1,1,1];
    col1 = [linspace(neg(1),mid(1),N1)',linspace(neg(2),mid(2),N1)',linspace(neg(3),mid(3),N1)'];
    col2 = [linspace(mid(1),pos(1),N2)',linspace(mid(2),pos(2),N2)',linspace(mid(3),pos(3),N2)'];
    col_black = [col1;col2];

    % - PLOT BAROTROPIC STREAM FUNCTION


    if isfield(opts, 'bstream') || isfield(opts, 'everything')
        figure(1);

        % - INTEGRATIONS - --------------------------------------------------
        % Barotropic streamfunction;
        PSIB = bstream(u*udim,zw*hdim,[y;ymax]*r0dim);

        img  = PSIB(2:end,:)';
        imgp = img; imgp(imgp<0)=NaN;
        imgn = img; imgn(imgn>-0)=NaN;

        contourf(RtD*x,RtD*(y),imgp,15,'LineStyle','none'); hold on;
        contourf(RtD*x,RtD*(y),imgn,15,'LineStyle','none'); hold on;

        plot_mask(summask,x,y); hold on

        contour(RtD*x,RtD*(y),imgp,15,'k'); hold on;
        contour(RtD*x,RtD*(y),imgn,15,'k--'); hold off;

        if fix_caxis
            caxis([opts.caxis_min,opts.caxis_max])
        else
            crange = max(abs(min(caxis)),abs(max(caxis)));
            caxis([-crange, crange]);
        end

        colorbar
        colormap(col_black)

        if plot_title
            title(['Barotropic Streamfunction (Sv) ', opts.title_add]);
        end

        xlabel('Longitude')
        ylabel('Latitude');


        if export_to_file
            exportfig(['bstream',opts.fname_add,'.eps'],10,[19,11])
        end

    end

    if isfield(opts, 'mstream') || isfield(opts, 'everything')
        figure(2);
        % PLOT OVERTURNING STREAMFUNCTION

        % Compute overturning streamfunction
        PSIG = mstream(v*udim,[x;xmax]*cos(yv(2:m+1))'*r0dim,zw*hdim);
        PSIG = [zeros(m+1,1) PSIG];

        PSIGp = PSIG; PSIGp(PSIGp<0)  = NaN;
        PSIGn = PSIG; PSIGn(PSIGn>0)  = NaN;
        contourf(RtD*([y;ymax+dy/2]-dy/2),zw*hdim',PSIGp',15); hold on
        contourf(RtD*([y;ymax+dy/2]-dy/2),zw*hdim',PSIGn',15,'--'); hold off
        colorbar
        cmin = min(min(PSIG(:,1:9)));
        cmax = max(max(PSIG(:,1:9)));
        fprintf('MOC+ = %f MOC- = %f MOC+ + MOC- = %f \n', max(PSIG(:)), ...
                min(PSIG(:)), max(PSIG(:)) + min(PSIG(:)))
        if plot_title
            title(['MOC (Sv) ',opts.title_add])
        end

        xlabel('latitude')
        ylabel('depth (m)')

        if fix_caxis
            caxis([opts.caxis_min,opts.caxis_max])
        else
            crange = max(abs(min(caxis)),abs(max(caxis)));
            caxis([-crange, crange]);
        end

        colormap(col_white)

        if export_to_file
            exportfig(['mstream',opts.fname_add,'.eps'],10,[19,10])
        end
    end

    if isfield(opts, 'everything')
        figure(3);
        % - CHECK SALINITY - ------------------------------------------------
        check = checksal(S,x,y,dfzt);
        vol   = sum(sum(1-surfm).*cos(y'))*dx*dy;
        fprintf(1,'Average salinity deficiency of %12.8f psu.\n', -check/vol)

        % Create Temperature and salinity
        % build longitudinal average over non-land cells
        Tl = zeros(m,l);
        Sl = zeros(m,l);

        for k = 1:l
            for j = 1:m
                count = 0;
                for i=1:n
                    if landm_int(i,j,k) == 0
                        count = count + 1;
                        Tl(j,k) = Tl(j,k) + T(i,j,k);
                        Sl(j,k) = Sl(j,k) + S(i,j,k);
                    else
                        S(i,j,k) = 0;
                        T(i,j,k) = 0;
                    end
                end
                Tl(j,k) = Tl(j,k) / count;
                Sl(j,k) = Sl(j,k) / count;
            end
        end

        % -------------------------------------------------------

        contourf(RtD*yv(1:end-1),z*hdim,Tl'+T0,15);
        %imagesc(RtD*yv(1:end-1),z*hdim,Tp'+T0);
        %set(gca,'ydir','normal');
        %pcolor(RtD*yv(1:end-1),z*hdim,Tp'+T0);
        colorbar
        title('Temperature')
        xlabel('Latitude')
        ylabel('z (m)')
        colormap(col_white)
        crange = max(abs(caxis))-T0;
        caxis([T0-crange, T0+crange]);

        if export_to_file
            exportfig('isothermals.eps',10,[20,7])
        end


        % -------------------------------------------------------
        figure(4);
        Tsurf = T(:,:,l);
        minT = T0+min(min(Tsurf));
        maxT = T0+max(max(Tsurf));

        img  = T0 + Tsurf';
        contourf(RtD*x,RtD*(y),img,20,'Visible', 'off'); hold on;
        set(gca,'color',[0.65,0.65,0.65]);

        contours = linspace(minT,maxT,20);
        imagesc(RtD*x,RtD*(y),img);

        hold off

        colorbar
        title('Surface Temperature', 'interpreter', 'none');
        xlabel('Longitude');
        ylabel('Latitude');
        colormap(col_white)

        crange = max(abs(min(caxis)),abs(max(caxis)))-T0;
        caxis(T0+[-crange, crange]);


        if export_to_file
            exportfig('sst.eps',10,[50,25])
        end

        figure(5);


        contourf(RtD*yv(1:end-1),z*hdim,Sl'+S0,15);
        %imagesc(Sp'+S0);
        %set(gca,'ydir','normal')
        %pcolor(RtD*yv(1:end-1),z*hdim,Sp'+S0);
        colorbar
        title('Isohalines')
        xlabel('Latitude')
        ylabel('z (m)')

        colormap(col_white)
        crange = max(abs(min(caxis)),abs(max(caxis)))-S0;
        caxis(S0+[-crange, crange])


        if export_to_file
            exportfig('isohalines.eps',10,[20,7])
        end

    end

end
