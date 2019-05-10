function [titles, cdata,h,e] = plot_cdata(fname, opts)
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
    
    if isfield(opts,'reduce_size')
        reduce_size = opts.reduce_size;
    else
        reduce_size = false;
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
    
    if plot_entry ~= 0 
        assert(plot_entry >= 2)
        assert(plot_entry <= size(cdata,2))
    end
    
    for i = 2:size(cdata,2)

        if plot_all || i == plot_entry

            if reduce_size
                
                dx = cdata(2:end,1)-cdata(1:end-1,1);
                
                sgn   = sign(dx);
                folds = sgn(2:end) ~= sgn(1:end-1);
                
                dy  = cdata(2:end,i)-cdata(1:end-1,i);
                dd  = sqrt(dx.^2+dy.^2);
                th  = max(dd)/0.8;
                row = 1;
                
                si = sign(dy(1)/dx(1));
                
                while row < size(cdata,1)-1
                    dx  = cdata(row+1,1)-cdata(row,1);
                    dy  = cdata(row+1,i)-cdata(row,i);
                    dd  = sqrt(dx.^2+dy.^2);

                    if (dd < th) && (~folds(row))
                        cdata(row+1,:) = [];
                        folds(row,:) = [];
                    else
                        row = row + 1;
                    end

                end
                
                % if rows > mx
                %     cdata = cdata(1:floor(rows/mx):rows,:);
                % end
                                 
            end
            
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
                e  = plot(xpoint,ypoint,'o','markerfacecolor',[1,.2,.2],'color',[1,1,1],'markersize',8); 
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