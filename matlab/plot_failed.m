function [files] = plot_failed(base_names, level, fig_ctr)

    comparefiles = false;
    if ~iscell(base_names)
        base_names   = {base_names};
    end
    
    if iscell(base_names) && numel(base_names) == 2
        comparefiles = true;
    end
    
    if comparefiles
        [out, list1] = system(['ls ', base_names{1}, '*']);
        [out, list2] = system(['ls ', base_names{2}, '*']);
    else
        [out, list1] = system(['ls ', base_names{1}, '*']);
    end

    [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = ...
        readfort44('fort.44');

    RtD   = 180/pi;   %[-] Radians to degrees

    if nargin < 2
        level = l;
    end

    if nargin < 3
        fig_ctr = 1;
    end

    % Collect files
    spaces1 = strfind(list1, ' ');
    nfiles1 = numel(spaces1) / 2 + 1;

    if comparefiles
        spaces2 = strfind(list2, ' ');
        nfiles2 = numel(spaces2) / 2 + 1;
        
        if nfiles1 ~= nfiles2
            fprintf('unequal amount of files in comparison\n')
            return;
        end       
    end
    
    files1  = cell(3,1);
    files2  = cell(3,1);
    sinds   = [1,spaces1(2:2:end)+1];
    einds   = [spaces1(1:2:end)-1,numel(list1)-1];

    for i = 1:nfiles1
        files1{i} = list1(sinds(i):einds(i));
        if comparefiles
            files2{i} = list2(sinds(i):einds(i));            
        end            
    end

    % assuming standard ordering
    NUNS = [6,3,4];
    LS   = [l,1,1];
    LVLS = [level,1,1];
    XTOT = [];
    Xnrm = [];
    
    fprintf('\n');
    for i = 1:nfiles1
        X    = readhdf5(files1{i});
        if comparefiles
            X = X - readhdf5(files2{i});
        end
        
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