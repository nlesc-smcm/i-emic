function [files] = plot_failed(base_name, level, fig_ctr)
    [out, list] = system(['ls ', base_name, '*']);

    [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = ...
        readfort44('fort.44');

    RtD   = 180/pi;              %[-]     Radians to degrees

    if nargin < 2
        level = l;
    end

    if nargin < 3
        fig_ctr = 1;
    end

    % Collect files
    spaces = strfind(list, ' ');
    nfiles = numel(spaces) / 2 + 1;
    files  = cell(3,1);
    sinds  = [1,spaces(2:2:end)+1];
    einds  = [spaces(1:2:end)-1,numel(list)-1];

    for i = 1:nfiles
        files{i} = list(sinds(i):einds(i));
    end

    % assuming standard ordering
    NUNS = [6,3,4];
    LS   = [l,1,1];
    LVLS = [level,1,1];
    XTOT = [];
    Xnrm = [];
    
    fprintf('\n');
    for i = 1:nfiles
        X    = readhdf5(files{i});
        XTOT = [XTOT;X];
        fprintf('-------------------------\n');      
        if i == 2
            P = X(end);
            X = X(1:end-1);
            fprintf(' P = %e\n', P);
        elseif i == 3
            G = X(end);
            X = X(1:end-1);
            fprintf(' G = %e\n', G);
        end
        fprintf('\n ||X%d|| = %e: \n', i, norm(X));
        X = reshape(X, NUNS(i),n,m,LS(i));                                              
        
        for xx = 1:NUNS(i)
                    
            img = squeeze(X(xx,:,:,LVLS(i)))';
            fprintf('     fig %d, ||xx||%d = %e \n', fig_ctr, xx, norm(img(:)));
            figure(fig_ctr); fig_ctr = fig_ctr+1;
            imagesc(RtD*x,RtD*(y), img);
            title(['model ', num2str(i), ', xx = ', num2str(xx), ...
                   ', level = ', num2str(LVLS(i))]);
            
            set(gca,'ydir','normal');
            colorbar;
        end
        fprintf('\n');
    end
    fprintf('--------------------------------\n');
    fprintf('\n  ||X|| = %e\n\n',  norm(XTOT));
end