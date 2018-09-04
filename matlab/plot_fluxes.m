function [out] = plot_fluxes(struct, fig_ctr, titlepre, opts)
    
    if nargin < 3
        titlepre = ''
    end
    
    if nargin < 4
        opts = [];
    end
    
    if isfield(opts, 'invert')
        invert = opts.invert;
    else
        invert = false;
    end      
    
    % find flux fields in struct
    RtD = 180/pi;
    [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');

    names = fieldnames(struct);
    ctr = 0;
    out = [];
    for i = 1:numel(names);
        match = regexp(names(i), '.*Flux', 'match');
        if numel(match{1}) > 0
            name   = match{1}{1};
            values = getfield(struct, name);
            if norm(values) == 0
                continue;
            end
            values = reshape(values, n, m);
            out = setfield(out, name, values);
            fprintf('%s\n', name)
            figure(fig_ctr); fig_ctr = fig_ctr + 1;
            imagesc(RtD*x,RtD*(y), values');
            set(gca, 'ydir', 'normal')
            cmap = my_colmap(caxis,0);
            colormap(cmap)
            colorbar
            xlabel('Longitude')
            ylabel('Latitude')
            title([titlepre, ' ', name])
            exportfig([name,'.eps'],10,[14,10],invert);
        end
    end
end