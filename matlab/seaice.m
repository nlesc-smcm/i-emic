function [] = seaice()
    
    global xmin xmax ymin ymax RtD n m nun x y dx dy 
    global t0 tvar s0 svar q0 qvar ta0 tavar 
    global sst, sss, tatm, qatm
    
    % specify phsyicsl dimensions
    xmin = 286 * pi / 180;
    xmax = 350 * pi / 180;
    ymin =  10 * pi / 180;
    ymax =  80 * pi / 180;

    % specify grid size
    n = 32;
    m = 16;
    
    % number of unknowns
    nun = 4;
    dim = nun*n*m;
    
    % constants
    t0   = 15;
    tvar = 10;
    s0   = 35;
    svar = 1;
    q0   = 8e-3;
    qvar = 3e-3;
    RtD  = 180/pi;
    
    % create grid
    grid();    
    
    sst  = idealizedTemp(t0, tvar);    
    sss  = idealizedSalt(s0, svar);
    tatm = idealizedTemp(t0-1, tvar);
    qatm = idealizedTemp(q0, qvar);
    
    x = zeros(dim, 1);
    
    F = rhs(x);   

end

function [F] = rhs(x)
    
    [H, Qtsa, Tsi, Msi] = extractsol(x);
    
    for j = 1:m
        for i = 1:n
            for XX = 1:nun
                row = find_row(i,j,XX);
                switch XX
                  case 1
                    val = (zeta * Tf(sss(i,j)) - sst - Qtsa(i,j)) / ...
                          ( rhoi * Lf ) - ... 
                          ( rhoo / rhoi) * ... 
                          (E0 + dEdT * Tsi(i,j) + dEdq * qatm(i, ...
                                                              j));
                  case 2
                    val = ..........%%%%%todotodotodo
                end
                
                F(row) = val
            end
        end
    end    
end


function [row] = find_row(i,j,XX)
    global n m nun
    
    row = nun * ( (j-1) * n  + (i-1) ) + XX;
end

function [H, Qtsa, Tsi, Msi] = extractsol(x)

    global n m nun 
    
    H    = zeros(n,m);
    Qtsa = zeros(n,m);
    Tsi  = zeros(n,m);
    Msi  = zeros(n,m);
    
    H(:)    = x(1:nun:end);
    Qtsa(:) = x(2:nun:end);
    Tsi(:)  = x(3:nun:end);
    Msi(:)  = x(4:nun:end);    
end      

function [] = grid()

    global xmin xmax ymin ymax n m x y dx dy
    
    dx = (xmax-xmin)/n;
    dy = (ymax-ymin)/m;
    
    x = ((1:n)-0.5)*dx + xmin;
    y = ((1:m)-0.5)*dy + ymin;
    
end

function [out] = idealizedTemp(mn, var)
    
    global xmin xmax ymin ymax n m x y 
    
    out = zeros(n,m);
    out = mn + var * repmat(cos(pi*y/ymax), [n, 1]);
end

function [out] = idealizedSalt(mn, var)
    
    global xmin xmax ymin ymax n m x y 
    
    out = zeros(n,m);
    out = mn + var * repmat(cos(pi*y/ymax) ./ cos(y), [n, 1]);
end


function [] = myvis(field, ttl)

    global x y RtD
    
    if nargin < 2 
        ttl = '';
    end
    
    imagesc(RtD*x,RtD*(y),field');
    set(gca,'ydir','normal');
    colorbar;
    col = my_colmap([min(field(:)), max(field(:))]);
    colormap(col);
    title(ttl)
end