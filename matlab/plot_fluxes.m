function [out] = plot_fluxes(struct, fig_ctr, titlepre)

    if nargin < 3
        titlepre = ''
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
            figure(fig_ctr) fig_ctr = fig_ctr + 1;
            set(gca, 'ydir', 'normal')
            cmap = my_colmap(caxis,0);
            colormap(cmap)
            colorbar
            xlabel('Longitude')
            ylabel('Latitude')
            title([titlepre, ' ', name])
        end
    end
end