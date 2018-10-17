function [titles, cdata] = plot_cdata(fname, opts)
    plot_fancy=false;
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
                cmax = max(cdata(:,i));
                cmn  = mean(cdata(:,i));
                cmap = jet(cmax);
                %cmap = my_colmap([0,cmax], -1, cmax);
                clrs = cmap(cdata(:,i),:);
                break;            
            end
        end
    end
    
    
    
    for i = 2:size(cdata,2)
        figure(i)

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
            uistack(s,'down');
        end
        title(titles{i})
        grid on;
    end
    
    if ((maxPsi > 0) && (minPsi > 0))
        figure(i+1);
        plot(cdata(:,1),cdata(:,maxPsi)+cdata(:,minPsi),lsty);
        title('max(Psi)+min(Psi)');
    end       

end