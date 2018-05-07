function [nrm] = plot_failed_residual(fname, maskname, level)

    if nargin < 1
        fname = 'failed_rhs.h5';
    end
    if nargin < 2
        maskname = 'fort.44';
    end

    if isempty(maskname)
        
    else

    end
    
    [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = ...
        readfort44(maskname);
    
    land = landm(2:end-1,2:end-1,2:l+1);

    if nargin < 3
        level = l;
    end

    fprintf(' using %s and %s at level %d\n', fname, maskname, level);

    RtD   = 180/pi;              %[-]     Radians to degrees

    res  = readhdf5(fname);
    size(res)
    
    oceanRes = true;
    atmosRes = false;    
    
    if (numel(fname) > 8)
        if strcmp(fname(end-7:end), 'first.h5')
            fprintf(' ocean residual\n');
            oceanRes = true;
            atmosRes = false;
        elseif strcmp(fname(end-8:end), 'second.h5')
            l   = 1;
            nun = 2;
            fprintf(' atmos residual\n');
            oceanRes = false;
            atmosRes = true;        
        end
    end

    sol = zeros(nun,n,m,l+la);
    idx = 1;
    for k = 1:l+la
        for j = 1:m
            for i = 1:n
                for XX = 1:nun
                    sol(XX,i,j,k) = res(idx);
                    idx = idx + 1;
                end
            end
        end
    end


    if (oceanRes)
        [u,v,w,p,T,S] = extractsol(sol);

        figure(1);
        mx = max(max(abs(u(:,:,level))'));
        img = -mx*2*land(:,:,level)' + u(:,:,level)';
        imagesc(RtD*x,RtD*(y), img);
        set(gca,'ydir','normal');

        colormap(parula)
        xlabel('Longitude')
        ylabel('Latitude');
        colorbar
        title('u')

        figure(2);
        mx = max(max(abs(v(:,:,level))'));
        img = -mx*2*land(:,:,level)' + v(:,:,level)';
        imagesc(RtD*x,RtD*(y), img);
        colorbar
        title('v')
        set(gca,'ydir','normal');
        colormap(parula)
        xlabel('Longitude')
        ylabel('Latitude');

        figure(3);
        mx = max(max(abs(w(:,:,level))'));
        img = -mx*2*land(:,:,level)' + w(:,:,level)';
        imagesc(RtD*x,RtD*(y), img);
        colorbar
        title('w')
        set(gca,'ydir','normal');
        colormap(parula)
        xlabel('Longitude')
        ylabel('Latitude');

        figure(4);
        mx  = max(max(abs(p(:,:,level))'));
        img = -mx*2*land(:,:,level)' + p(:,:,level)';
        imagesc(RtD*x,RtD*(y), img);
        colorbar
        title('p')
        set(gca,'ydir','normal');
        colormap(parula)
        xlabel('Longitude')
        ylabel('Latitude');

        figure(5);
        mx  = max(max(abs(T(:,:,level))'));
        img = -mx*2*land(:,:,level)' + T(:,:,level)';
        imagesc(RtD*x,RtD*(y), img);
        colorbar
        title('T')
        set(gca,'ydir','normal');
        colormap(parula)
        xlabel('Longitude')
        ylabel('Latitude');
         
        figure(6);
        mx  =  max(max(abs(S(:,:,level))'));

        %S(~logical(land)) = S(~logical(land)) - .2;
        
        img = S(:,:,level)';
        img = img - 2*max(abs(img(:)))*land(:,:,level)';

        % img(img==0)=NaN;
        imagesc(RtD*x,RtD*(y), img);
        colorbar
        title('S')
        set(gca,'ydir','normal');

        xlabel('Longitude')
        ylabel('Latitude');

        fprintf('|u| = %e\n',norm(u(:)));
        fprintf('|v| = %e\n',norm(v(:)));
        fprintf('|w| = %e\n',norm(w(:)));
        fprintf('|p| = %e\n',norm(p(:)));
        fprintf('|T| = %e\n',norm(T(:)));
        fprintf('|S| = %e\n',norm(S(:)));
        fprintf('---------\n');
    
    elseif atmosRes
        T = squeeze(sol(1,:,:,:));
        q = squeeze(sol(2,:,:,:));
        P = res(n*m*nun+1)
        
        figure(7) 
        mx  =  max(max(abs(T(:,:,1))'));
        img = -mx*2*land(:,:,size(land,3))' + T(:,:,1)';
        imagesc(RtD*x,RtD*(y), img);
        colorbar
        title('T')
        set(gca,'ydir','normal');
        colormap(parula)
        xlabel('Longitude')
        ylabel('Latitude');

        figure(8) 
        mx  =  max(max(abs(q(:,:,1))'));
        img = -mx*2*land(:,:,size(land,3))' + q(:,:,1)';
        imagesc(RtD*x,RtD*(y), img);
        colorbar
        title('q')
        set(gca,'ydir','normal');
        colormap(parula)
        xlabel('Longitude')
        ylabel('Latitude');


        fprintf('|T| = %e\n',norm(T(:)));
        fprintf('|q| = %e\n',norm(q(:)));
        fprintf('|P| = %e\n',norm(P(:)));
        
    end
    
    nrm = norm(res(:));
    fprintf('total |X| = %e\n', nrm);
    
end
