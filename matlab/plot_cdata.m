function [titles, cdata,h] = plot_cdata(fname, opts)
    if nargin < 2 
        opts = [];
    end
    
    if isfield(opts,'colors')
        plot_fancy=true;
        coldata = opts.colors;
    else
        plot_fancy=false;
    end
    
    if isfield(opts,'lsty')
        lsty = opts.lsty;        
    else
        lsty={'k.-'};
    end

    if isfield(opts,'prepend') && isfield(opts,'prefiles')
        npfiles  = opts.prepend;
        prefiles = opts.prefiles;
    else
        npfiles = 0;
    end
    
    if isfield(opts,'hold')
        holdfig = opts.hold;
    else
        holdfig = false;
    end

    if isfield(opts,'exportnum')
        exportnum = opts.exportnum;
    else
        exportnum = 0;
    end

    if isfield(opts,'plot_all')
        plot_all = opts.plot_all;
    else
        plot_all = true;
    end

    if isfield(opts,'plot_entry')
        plot_entry = opts.plot_entry;
        plot_all = false;
    else
        plot_entry = 0;
        plot_all = true;
    end

    if isfield(opts,'invert')
        invert = opts.invert;
    else
        invert = false;
    end

    if isfield(opts,'partrans')
        partrans = opts.partrans;
    else
        partrans = @(x) x;
    end

    if isfield(opts,'parname')
        parname = opts.parname;
    else
        parname = 'par';
    end
    
    if isfield(opts,'point')
        point = opts.point;
    else
        point = Inf;
    end

    if isfield(opts,'xzmpoint') && isfield(opts,'yzmpoint')
        xzmpoint = opts.xzmpoint;
        yzmpoint = opts.yzmpoint; 
    else
        xzmpoint = 0;
        yzmpoint = 0;
    end
    
    if (nargin < 1)
        fname = 'cdata.txt'
    end

    cdata = [];
    for i = 1:npfiles
        [t,c] = load_cdata(prefiles{i});
        cdata = [cdata; c];
    end       
    
    [titles, c] = load_cdata(fname);
    cdata = [cdata; c];

    maxPsi = 0;
    minPsi = 0;
    
    if plot_fancy
        for i = 1:size(cdata,2)
            if strcmp(titles{i}, coldata)
                dmax = max(abs(cdata(:,i)));
                cmap = jet(129);
                clrs = cmap(ceil(128*abs(cdata(:,i))/dmax),:);
                break;            
            end
        end
    end    
    
    for i = 2:size(cdata,2)

        if plot_all || i == plot_entry
            
            figure(i)
            if holdfig
                hold on;
            end

            if strcmp(titles{i}, 'max(Psi)')
                maxPsi = i;
            elseif strcmp(titles{i}, 'min(Psi)')
                minPsi = i;
            end
            
            if strcmp(titles{i}, '||F||') || strcmp(titles{i}, 'ds') || ...
                    strcmp(titles{i}, 'Tol')
                cdata(:,i) = abs(cdata(:,i));
                h = semilogy(partrans(cdata(:,1)),cdata(:,i),lsty{:});
                xlabel(parname);
            else
                h = plot(partrans(cdata(:,1)),cdata(:,i),lsty{:});
                xlabel(parname);
            end
            if plot_fancy
                hold on
                s = scatter(partrans(cdata(:,1)),cdata(:,i), 2, clrs, 'o','filled'); 
                uistack(s,'bottom');
                hold off;
            end

            if point > 0
                if point > size(cdata,1)
                    point = size(cdata,1);
                end
                xpoint = partrans(cdata(point,1));
                ypoint = cdata(point,i);
                hold on
                e  = plot(xpoint,ypoint,'o','markerfacecolor','w'); 
                hold off;
            end
            
            if (xzmpoint >  0) && (yzmpoint >  0)
                xlim([xpoint-xzmpoint,xpoint+xzmpoint])
                ylim([ypoint-yzmpoint,ypoint+yzmpoint])
            end
            
            title(titles{i})
            grid on;
            
            if holdfig
                hold off;
            end
            
            if exportnum==0 && invert
                invertcolors();
            end

            if i == exportnum
                %legend([s,e],coldata,'location','northwest');                
                exportfig([titles{i},'.eps'], 9, [10,10], invert);
            end
        end
    end
    
    if ((maxPsi > 0) && (minPsi > 0))    
        figure(i+1);
        
        if holdfig
            hold on;
        end
        
        plot(partrans(cdata(:,1)),cdata(:,maxPsi)+cdata(:,minPsi),lsty{:});
        title('max(Psi)+min(Psi)');
        
        if holdfig
            hold off;
        end
    end       
    
end