labels = {};
fnames = {};
dnames = {};

% create avi file
createAvi = false;

% add interpolated frames
addInterp = false;

Ntopos = 14;
for i = 1:Ntopos
    labels = [labels, [num2str(70-5*i) 'Ma']];
    fnames = [fnames, ['state_topo_', num2str((i-1)), '.h5']];
    dnames = [dnames, ['mask_', num2str((i-1)), '.mask']];
end
fnames{1} = 'ocean_input.h5';

fname = 'paleomovie_interp.avi';
if createAvi
    writerObj = VideoWriter(fname, 'Motion JPEG AVI');
    writerObj.FrameRate = 10;
    writerObj.Quality = 90;
    open(writerObj);
    fhandle = figure('units','pixels','position',[100,100,1000,600]);
    set(gca,'position',[0.05 0.1 .92 0.85],'units','normalized');
    set(gca,'color','w','fontsize',15);
end

counter = 0;
clear opts;
clear opts_interp;

for i = 1:Ntopos
    counter = counter + 1;

    opts.bstream   = true;
    if i <= Ntopos
        opts.title_add = [labels{i}];
    end
    opts.fname_add = ['-inv-',num2str(counter)]

    opts.fix_caxis = true;
    opts.caxis_min = -40;
    opts.caxis_max = 40;
    opts.invert = true;

    fprintf('\n   %d %d --- \n', counter, i)

    plot_ocean(fnames{i}, dnames{i}, opts)
    
    if createAvi
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
    end

    r = sin(pi*linspace(0,1,30)/2).^2;
    r = linspace(0,1,20);
    
    if i < Ntopos && addInterp
        for k = r(2:end-1)
            fprintf('\n     --    %f \n', k)

            counter = counter + 1;
            opts_interp.bstream    = true;
            opts_interp.title_add  = [labels{i},'->',labels{i+1}];
            opts_interp.fname_add  = ['-inv-',num2str(counter)];
            opts_interp.solfile2   = fnames{i+1};
            opts_interp.maskfile2  = dnames{i+1};
            opts_interp.interp_par = k;
            opts_interp.fix_caxis = true;
            opts_interp.caxis_min = -40;
            opts_interp.caxis_max = 40;


            plot_ocean(fnames{i}, dnames{i+1}, opts_interp)
            if createAvi
                frame = getframe(gcf);
                writeVideo(writerObj, frame);
            end
        end
    end
end
if createAvi
    close(writerObj);
end