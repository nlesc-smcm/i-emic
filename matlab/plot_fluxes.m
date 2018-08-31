function [] = plot_fluxes(struct, figid)
% find flux fields in struct
    RtD = 180/pi;    
    [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');
        
    names = fieldnames(struct);
    ctr = 0;
    for i = 1:numel(names);
        match = regexp(names(i), '.*Flux', 'match');
        if numel(match{1}) > 0
            name   = match{1}{1};
            values = getfield(struct, name);
            values = reshape(values, n, m);
            figure(figid+ctr)
            imagesc(RtD*x,RtD*(y), values');
            set(gca, 'ydir', 'normal')
            cmap = my_colmap(caxis,0);
            colormap(cmap)
            colorbar
            xlabel('Longitude')
            ylabel('Latitude')
            title(name)
            ctr = ctr + 1;
        end
    end

end