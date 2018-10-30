function [titles, cdata] = plot_cdata(fname, opts)
    plot_fancy = false;
    holdfig = false;
    npfiles=0;
    lsty = 'k.-';
    
    if nargin > 1
        if isfield(opts,'colors')
            plot_fancy=true;
            coldata = opts.colors;
        end
        
        if isfield(opts,'lsty')
            lsty = opts.lsty;
        end

        if isfield(opts,'prepend') && isfield(opts,'prefiles')
            npfiles  = opts.prepend;
            prefiles = opts.prefiles;
        end
        
        if isfield(opts,'hold')
            holdfig = opts.hold;
        end
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
            semilogy(cdata(:,1),cdata(:,i),lsty);
            xlabel('par');
        else
            plot(cdata(:,1),cdata(:,i),lsty);
            xlabel('par');
        end
        if plot_fancy
            hold on
            s = scatter(cdata(:,1),cdata(:,i), 20, clrs, 'o','filled'); hold off
            uistack(s,'bottom');
        end
        title(titles{i})
        grid on;
        
        if holdfig
            hold off;
        end
    end
    
    if ((maxPsi > 0) && (minPsi > 0))
        figure(i+1);

        if holdfig
            hold on;
        end

        plot(cdata(:,1),cdata(:,maxPsi)+cdata(:,minPsi),lsty);
        title('max(Psi)+min(Psi)');
        
        if holdfig
            hold off;
        end
    end       

end