function [] = seaice()
    
    global xmin xmax ymin ymax RtD n m nun x y dx dy 
    global t0 tvar s0 svar q0 qvar ta0 tavar 
    global sst, sss, tatm, qatm
    
    % specify physical dimensions
    RtD  = 180 / pi;
    xmin = 286 / RtD;
    xmax = 350 / RtD;
    ymin =  10 / RtD;
    ymax =  80 / RtD;

    % specify grid size
    n = 32;
    m = 16;
    
    % number of unknowns
    nun = 4;
    dim = nun*n*m;
    
    % general physical constants
    t0o  = 15;
    tvar = 10;
    t0i  = -15;
    s0   = 35;
    svar = 1;
    q0   = 8e-3;
    qvar = 3e-3;
    
    % ice formation constants
    ch   = 0.0058;     % empirical constant
    utau = 0.02;       % ms^{-1}, skin friction velocity
    rhoo = 1.024e3;    % kg m^{-3}, sea water density
    rhoi = 0.913e3;    % kg m^{-3}, ice density
    cpo  = 4.2e3;      % J kg^{-1} K^{-1], sea water heat capacity
    Lf   = 3.34e5;     % J kg^{-1}, latent heat of fusion of ice 
    
    zeta = ch * utau * rhoo * cpo;
    
    % sublimation constants
    % Parameters for saturation humidity over ocean and ice
    c1 = 3.8e-3;  % (kg / kg)
    c2 = 21.87;   %
    c3 = 265.5;   % (K)
    c4 = 17.67;   %
    c5 = 243.5;   % (K)
    
    % Calculate background saturation specific humidity according to
    % [Bolton,1980], T in \deg C
    qsi = c1 * exp(c2 * t0i / (t0i + c3));
    
    % create grid
    grid();    
    
    sst  = idealizedTemp(t0o, tvar);    
    sss  = idealizedSalt(s0, svar);
    tatm = idealizedTemp(t0o-1, tvar);
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
                          (E0 + dEdT * Tsi(i,j) + dEdq * qatm(i,j));
                  case 2
                    val = -Qtsa(i,j) + ...
                          (sigma0 / 4) * S(y(j)) * (1-alpha) * C0 - ...
                          mu * (Tsi(i,j) - tatm(i,j)) - ...
                          rhoo * Ls * (E0 + dEdT * Tsi(i,j) + dEdq * qatm(i,j));
                  case 3 
                    val = Tsi(i,j) - Tf(sss(i,j)) - ...
                          Qtsa(i,j) * H(i,j) / Ic;
                  case 4
                    val = Msi(i,j) - ...
                          (1/2) * (1 + tanh((H(i,j) - taus) * epsilon));                    
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