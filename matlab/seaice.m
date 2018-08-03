function [X, J, F] = seaice()

    global xmin xmax ymin ymax RtD n m nun aux x y dx dy
    global IC
    global t0o t0a t0i Q0 Qvar H0 M0
    global tvar s0 svar q0 qvar ta0 tavar

    % forcing from external models
    global sst sss qatm patm tatm

    % model parameters and functions
    global zeta a0 Tf Lf rhoi rhoo E0i E0o dEdT dEdq
    global sun0 S Ls alpha c0 muoa
    global Ic
    global taus epsilon
    global Ti
    global nus qdim tdim dqso dqsi eta

    global combf solf maskf latf % continuation parameters

    combf = 1.0;
    solf  = 1.0;
    maskf = 1.0;
    latf  = 1.0;

    % specify physical dimensions
    RtD  = 180 / pi;
    xmin = 286 / RtD;
    xmax = 350 / RtD;
    ymin =  10 / RtD;
    ymax =  80 / RtD;

    % specify grid size (2deg)
    n = 8;
    m = 8;

    % create grid
    grid();

    % number of unknowns
    aux = 1;
    nun = 4;
    dim = nun*n*m;
    len = dim + aux;

    taus    = 0.01;     % m, threshold ice thickness
    epsilon = 1e-2;     % Heavyside approximation steepness

    % sublimation constants, parameters for saturation humidity over ocean
    % and ice
    c1 = 3.8e-3;   % (kg / kg)
    c2 = 21.87;    %
    c3 = 265.5;    % (K)
    c4 = 17.67;    % (K)
    c5 = 243.5;    % (K)

    % general physical constants
    t0o  =  15;
    t0a  =  15;
    t0i  =  c3*c4*t0o / (c2*c5+(c2-c4)*t0o);
    tvar =  15;
    s0   =  35;
    svar =  1;
    q0   =  8e-3;
    qvar =  1e-3;
    H0   =  taus;

    % ice formation constants
    ch   = 0.0058;     % empirical constant

    utau = 0.02;       % ms^{-1}, skin friction velocity
    rhoo = 1.024e3;    % kg m^{-3}, sea water density
    rhoi = 0.913e3;    % kg m^{-3}, ice density
    rhoa = 1.25;       % kg m^{-3}, atmospheric density
    cpo  = 4.2e3;      % W s kg^{-1} K^{-1], sea water heat capacity
    Lf   = 3.347e5;    % J kg^{-1}, latent heat of fusion of ice
    Ls   = 2.835e6;    % J kg^{-1}, latent heat of sublimation of ice
    Ic   = 2.166;      % W m^{-1} K^{-1}, constant ice conductivity

    % combined parameter
    zeta = ch * utau * rhoo * cpo

    ce  = 1.3e-03; % Dalton number
    uw  = 8.5;     % ms^{-1}, mean atmospheric surface wind speed

    eta = ( rhoa / rhoo ) * ce * uw;

    % I-EMIC nondim components
    qdim = qvar;
    tdim = 1.0;
    nus  = eta * qdim;

    % Calculate background saturation specific humidity according to
    % [Bolton,1980], T in \deg C
    qsi = @(t0i) c1 * exp(c2 * t0i ./ (t0i + c3));
    qso = @(t0o) c1 * exp(c4 * t0o ./ (t0o + c5));

    % Shortwave radiation constants and functions
    alpha = 0.3;      % albedo
    sun0  = 1360;     % solar constant
    c0    = 0.43;     % atmospheric absorption coefficient
    Ch    = 1.22e-3;  % a constant...
    cpa   = 1000;     % heat capacity

    % exchange coefficient
    muoa  = rhoa * Ch * cpa * uw;

    % latitudinal dependence shortwave radiation
    S = @(y) (1 - .482 * (3 * (sin(y)).^2 - 1.) / 2.);

    % freezing temperature (dominant term)
    a0 = -0.0575;
    Tf = @(S) a0 * (S + s0);

    % Background sublimation and derivatives
    E0i   =  eta * ( qsi(t0i) - q0 );

    E0o   =  eta * ( qso(t0o) - q0 );

    dqso  = (c1 * c4 * c5) / (t0o + c5).^2 * ...
            exp( (c4 * t0i) / (t0o + c5) );

    dqsi  = (c1 * c2 * c3) / (t0i + c3).^2 * ...
            exp( (c2 * t0i) / (t0i + c3) );

    dEdT  =  nus * tdim / qdim * dqsi;
    dEdq  =  nus * -1;

    % Background heat flux and variation
% $$$     Q0 = zeta*(Tf(0) - t0o) - rhoo * Lf * E0i;
% $$$
% $$$     Q0 = sun0 / 4 * S(mean(y)) * (1-alpha) * c0 - muoa * 2 - rhoo * ...
% $$$          Ls * E0i;
% $$$
% $$$     Q0 = sun0 / 4 * S(mean(y)) * (1-0.8) * c0 - muoa * 10 - rhoo * ...
% $$$          Ls * E0i;

    Q0 = -100;

    Qvar = zeta;

    % ice surface temperature (linearized)
    Ti = @(Q,H,S) (Tf(S) - t0i) + (Q0*H0 + H0*Qvar*Q ...
                                   + Q0*H) / Ic;

    % create integral coefficients
    IC = intcoeff();
    A  = sum(IC);
    fprintf(' total area = %2.12f\n', A);

    % initialize forcing and state
    sst  = idealizedTemp(0, tvar)-10;
    sss  = idealizedSalt(0, svar);
    tatm = idealizedTemp(0, tvar)-10;
    qatm = idealizedTemp(0, qvar);
    patm = ones(n,m);

    % creating test values
    fid = fopen('testResults', 'w');
    for par = [1.0, 0.1]
        combf = par;
        for scale = [0.0, 1.0, 1.234]
            X = scale*ones(len,1);
            F = rhs(X);
            fprintf('X = %1.3f, par = %1.1f, ||F||two = %1.12e\n', ...
                    scale, combf, norm(F));
            fprintf(fid, '%1.12e\n', norm(F));

        end
    end

    fprintf('\n');
    fprintf(fid, '\n');


    for par = [1.0, 0.1]
        combf = par;
        for scale = [0.0,1.0,1.234]
            X = scale*ones(len,1);
            J = jac(X);

            fprintf('X = %1.3f, par = %1.1f, ||J||inf = %1.12e\n', ...
                    scale, combf, norm(J,Inf));
            fprintf('X = %1.3f, par = %1.1f, ||J||one = %1.12e\n', ...
                    scale, combf, norm(J,1));
            fprintf('X = %1.3f, par = %1.1f, ||J||frb = %1.12e\n', ...
                    scale, combf, norm(J,'fro'));

            fprintf(fid, '%1.12e\n', norm(J,Inf));
            fprintf(fid, '%1.12e\n', norm(J,1));
            fprintf(fid, '%1.12e\n', norm(J,'fro'));
        end
    end
    fprintf('\n');
    return

    kmax = 10;
    ord = [];
    for i = 1:nun
        ord = [ord, i:nun:dim];
    end
    if (aux > 0)
        ord = [ord, ord(end)+1:len];
    end

    X = randn(len,1);
    F = rhs(X);
    J   = jac(X);
    Jn  = numjacob(@rhs, X);

% $$$     vsm(J(ord,ord));
% $$$     vsm(Jn(ord,ord));
% $$$     vsm(Jn(ord,ord)-J(ord,ord));


    FH = F(1:nun:dim);

    fprintf('   integral =  %e\n', IC'*FH );

    for i = 1:aux
        fprintf(' F(aux %d) = %e\n', i, F(dim+i))
    end

    for i = 1:kmax
        J   = jac(X);
        Jn  = numjacob(@rhs, X);

        dX = J \ -F;
        X  = X + dX;

        F  = rhs(X);
        fprintf('\n%2d| %e %e %e\n', i, norm(dX), norm(F), condest(J));
        if norm(dX) < 1e-12
            break;
        end

        [H,Q,M,T] = extractsol(X);
        FH = F(1:nun:dim);
        fprintf('   integral =  %e\n', IC'*(M(:).*FH) );

        for i = 1:aux
            fprintf('   F(aux %d) = %e\n', i, F(dim+i))
        end

    end

    [H,Q,M,T] = extractsol(X);
    [~, QSos, ~, EmiP] = fluxes(X);

% $$$     figure(1)
% $$$     Himg = (H'+H0);
% $$$     imagesc(RtD*x, RtD*y, Himg);
% $$$     set(gca,'ydir','normal'); colorbar
% $$$     title('H')

    figure(2)
    imagesc(RtD*x, RtD*y, M');
    set(gca,'ydir','normal'); colorbar
    title('M')

% $$$     figure(3)
% $$$     SST = sst'+t0o;
% $$$     imagesc(RtD*x, RtD*y, SST); set(gca,'ydir','normal'); colorbar
% $$$     title('SST')
% $$$     colormap(my_colmap(caxis, t0o));

% $$$     figure(4)
% $$$     Tsi = Tf(sss') + (Q0*H0+H0*Qvar*Q'+Q0*H') / Ic;
% $$$     imagesc(RtD*x,RtD*y,Tsi);
% $$$     set(gca,'ydir','normal');
% $$$     colorbar;
% $$$     title('Ti_2')

    figure(5)
    Tsi = T' + t0i;
    Tsi(M'<1e-2)=NaN;
    imagesc(RtD*x, RtD*y, Tsi);
    set(gca,'ydir','normal');
    colorbar;
    title('Ti')

% $$$     figure(6)
% $$$     Qimg = Qvar*Q' + Q0;
% $$$     imagesc(RtD*x,RtD*y,Qimg); set(gca,'ydir','normal'); colorbar;
% $$$     title('Q')
% $$$
% $$$     figure(7)
% $$$     imagesc(RtD*x,RtD*y,QSos'); set(gca,'ydir','normal'); colorbar;
% $$$     title('Q_S^{os}')
% $$$
% $$$     figure(8)
% $$$     imagesc(RtD*x,RtD*y,EmiP'); set(gca,'ydir','normal'); colorbar;
% $$$     title('E_S - P')


end

function [x] = initialsol()

    global m n nun aux
    dim = m*n*nun;
    x = zeros(dim,1);

    H = 1e-3 * randn(n,m);
    Q = 1e-3 * randn(n,m);
    T = 1e-3 * randn(n,m);

    M = zeros(n,m);
    M(end,:) = 1;

    x(1:nun:dim) = H(:);
    x(2:nun:dim) = Q(:);
    x(3:nun:dim) = M(:);
    x(4:nun:dim) = T(:);

    R = ones(aux,1);

    x = [x;R];

end

function [J,Al] = jac(x)

    global m n nun aux y

    global IC
    global sst sss qatm tatm patm

    global t0o t0a t0i Q0 Qvar H0

    global zeta Tf Lf rhoi rhoo E0i dEdT dEdq
    global sun0 S Ls alpha c0 muoa
    global Ic
    global taus epsilon
    global Ti
    global combf nus solf maskf latf % continuation parameter

    [H, Qtsa, Msi, ~, R] = extractsol(x);

    % define indices for unknowns
    HH = 1;
    QQ = 2;
    MM = 3;
    TT = 4;
    GG = 5;

    % initialize dependency grid
    Al = zeros(n, m, nun, nun+aux);

    % define dependencies

    % dHdt equation
    Al(:,:,HH,HH) = -rhoo * Lf / zeta * dEdT * Q0 / Ic;

    Al(:,:,HH,QQ) = -Qvar / zeta - ...
        rhoo * Lf / zeta * dEdT * H0 * Qvar / Ic;

    % Qtsa equation
    Al(:,:,QQ,QQ) = Qvar / muoa;
    Al(:,:,QQ,TT) = 1 + combf * latf * rhoo * Ls / muoa * dEdT;

    % Msi equation
    for j = 1:m
        for i = 1:n
            Al(i,j,MM,HH) = -(1 / 2) * ( ...
                tanh( (H(i,j) + 1e-6) / epsilon ) ...
                - tanh( H(i,j) / epsilon ) ) / 1e-6;
        end
    end

    Al(:,:,MM,MM) =  1;

    % Tsi equation
    Al(:,:,TT,HH) =  Q0 / Ic;
    Al(:,:,TT,QQ) =  H0 * Qvar / Ic;
    Al(:,:,TT,TT) = -1;

    [co, ico, jco] = assemble(Al);

    % auxiliary integrals
    if (aux > 0)

        [~,Qsos,~,EmiP] = fluxes(x);
        A  = sum(IC);

        for j = 1:m
            for i = 1:n
                sr = (j-1)*n+i;

                % dFGG / dQ
                row = find_row(i,j,GG);
                col = find_row(i,j,QQ);
                ico = [ico; row];
                jco = [jco; col];
                val = -1.0 * IC(sr) * Msi(sr) * Qvar / rhoo / Lf;
                co  = [co; val];

                % dFGG / dM
                row = find_row(i,j,GG);
                col = find_row(i,j,MM);
                ico = [ico; row];
                jco = [jco; col];
                val = IC(sr) * (Qsos(sr)-EmiP(sr));
                co  = [co; val];

                % dFGG / dT
                row = find_row(i,j,GG);
                col = find_row(i,j,TT);
                ico = [ico; row];
                jco = [jco; col];
                val = -1.0 * IC(sr)* Msi(sr) * dEdT;
                co  = [co; val];
            end
        end

        %dFG / dG
        row = find_row(i,j,GG);
        col = find_row(i,j,GG);

        ico = [ico; row];

        jco = [jco; col];
        co  = [co;  -A];

    end

    J = spconvert([ico,jco,co]);

end

function [F] = rhs(x)

    global m n nun y aux

    global IC

    global sst sss qatm tatm patm

    global t0o t0a t0i Q0 Qvar H0 M0

    global zeta Tf Lf rhoi rhoo E0i E0o dEdT dEdq
    global sun0 S Ls alpha c0 muoa
    global Ic
    global taus epsilon nus eta q0
    global Ti
    global combf solf maskf latf % continuation parameters

    [H, Qtsa, Msi, T, R] = extractsol(x);

    G = 1;

    N = size(x,1);
    F = zeros(N,1);

    for j = 1:m
        for i = 1:n
            for XX = 1:nun
                switch XX
                  case 1
                    
                    Tsi = Ti(Qtsa(i,j), H(i,j), sss(i,j) );

                    val =  Tf(sss(i,j)) - sst(i,j) - t0o - Q0 / zeta  ...
                           - Qvar / zeta * Qtsa(i,j) - ...
                           ( rhoo * Lf / zeta ) * ...
                           ( E0i + dEdT * Tsi + dEdq * ...
                             qatm(i,j) );

                  case 2
                    val = 1/muoa*Q0 + Qvar/muoa * Qtsa(i,j) - ...
                          (combf * solf * sun0 / 4 / muoa) * S(y(j)) * (1-alpha) * c0 + ...
                          (T(i,j) + t0i - tatm(i,j) - t0a ) + ...
                          (combf * latf * rhoo * Ls / muoa) * ...
                          (E0i + dEdT * T(i,j) + dEdq * ...
                           qatm(i,j));

                  case 3
                    val = Msi(i,j) - ...
                          (1 / 2) * (1 + tanh( (H(i,j) / ...
                                                            epsilon )));

                  case 4
                    val = Tsi - T(i,j);
                    
                end

                row = find_row(i,j,XX);
                F(row) = val;
            end
        end
    end

    if (aux > 0)
        [~, ~, ~, ~, R] = extractsol(x);
        [~,Qsos,~,EmiP,FluxDiff] = fluxes(x);

        A = sum(IC);
        IC'*FluxDiff(:);
        
        for aa = 1:aux
            switch aa
              case G
                value = IC' * (Msi(:) .* (Qsos(:) - EmiP(:))) ...
                        - R(G) * A;
            end
            F(row + aa) = value;
        end
    end
end

% single dof
function [QS, QSos, QSoa, EmiP, FluxDiff] = fluxes(x)

    global n m nun aux
    global sst sss qatm tatm patm
    global t0o t0a t0i Q0 Qvar H0 M0
    global zeta a0 Tf Lf rhoi rhoo E0i dEdT dEdq
    global sun0 S Ls alpha c0 muoa
    global Ic
    global taus epsilon
    global Ti
    global combf solf maskf   % continuation parameters
    global nus qdim tdim dqso dqsi eta

    len = n*m;
    QS  = zeros(n,m);

    [H, Qtsa, Msi, T] = extractsol(x);

    QSos     = zeros(n,m);
    QSoa     = zeros(n,m);
    EmiP     = zeros(n,m); % dimensional sublimation term
    FluxDiff = zeros(n,m); % dimensional sublimation term

    for j = 1:m
        for i = 1:n

            So = sss(i,j);
            To = sst(i,j);
            qs = Qtsa(i,j);
            qa = qatm(i,j);
            pa = patm(i,j);

            QSos(i,j) =  (                            ...  %
                zeta * combf * ( Tf(So) - (To+t0o) )  ...  % ! QTos component
                - ( Qvar * qs + combf * Q0 )) / rhoo / Lf; % ! QTsa component

            QSoa(i,j) = nus * ( ...
                (tdim / qdim) * dqso * To ...
                - qa - pa);

            EmiP(i,j) = nus * (tdim / qdim * dqsi * T(i,j) - qa - pa);

            % total salinity flux
            QS(i,j) = QSoa(i,j) + Msi(i,j) * (rhoo*Lf/zeta*QSos(i,j) ...
                                              - QSoa(i,j));

            FluxDiff(i,j) = Msi(i,j) * (QSos(i,j) - EmiP(i,j));
        end
    end
end


% integral coefficients for single dof surface row
function [IC] = intcoeff()
    global n m nun aux
    global dx dy x y
    IC = zeros(n,m);
    for j = 1:m
        for i = 1:n
            IC(i,j) = cos( y(j) ) * dx * dy;
        end
    end
    IC = IC(:);
end


function [co, ico, jco, beg] = assemble(Al)

    n    = size(Al,1);
    m    = size(Al,2);
    nunR = size(Al,3);
    nunC = size(Al,4);

    % This is an assemble where we assume all dependencies are
    % located at the same grid point. Otherwise we would have a
    % higher dimensional Al.

    co  = zeros(m*n*nunR*nunC,1); % values
    ico = zeros(m*n*nunR*nunC,1); % column indices
    jco = zeros(m*n*nunR*nunC,1); % column indices
    beg = zeros(m*n*nunR + 1, 1); % row pointer

    elm_ctr = 1; % element counter
    co_ctr  = 0;
    ico_ctr = 0;
    jco_ctr = 0;
    beg_ctr = 0;

    for j = 1:m
        for i = 1:n
            for A = 1:nunR

                % obtain row
                row = find_row(i,j,A);

                % fill beg with element counter
                beg_ctr = beg_ctr + 1;
                beg(beg_ctr) = elm_ctr;

                for B = 1:nunC
                    value = Al(i,j,A,B);

                    if (abs(value) > 0)

                        ico_ctr = ico_ctr + 1;
                        jco_ctr = jco_ctr + 1;
                        co_ctr  = co_ctr + 1;
                        elm_ctr = elm_ctr + 1;

                        % fill coefficient
                        co(co_ctr) = value;

                        % fill row index
                        ico(ico_ctr) = row;

                        % fill column index
                        col = find_row(i,j,B);
                        jco(jco_ctr) = col;
                    end
                end
            end
        end
    end

    % Truncate arrays
    co  = co(1:co_ctr);
    ico = ico(1:ico_ctr);
    jco = jco(1:jco_ctr);
    beg = beg(1:beg_ctr);
end

function [co, ico, jco] = assemble_fast(Al)

    global m n nun

    % This is an assemble where we assume all dependencies are
    % located at the same grid point. Otherwise we would have a
    % higher dimensional Al. The fast version does not compute beg.

    val = zeros(m*n*nun*nun,1);

    ctr = 0;
    for j = 1:m
        for i = 1:n
            for A = 1:nun
                for B = 1:nun
                    ctr = ctr + 1;
                    val(ctr) = Al(i,j,A,B);
                end
            end
        end
    end

    id  = logical(abs(val)>0);
    co  = val(id);

    ndim = m*n*nun;

    ico = repmat((1:ndim),[nun,1]);
    ico = ico(:);
    ico = ico(id);

    jco = (1:nun:ndim) + repmat((0:nun-1)',[nun,1]);
    jco = jco(:);
    jco = jco(id);
end

function [row] = find_row(i,j,XX)

    global n m nun aux

    dim = n*m*nun;
    len = dim + aux;

    if (XX > nun)
        row = dim + XX - nun;
    else
        row = nun * ( (j-1) * n  + (i-1) ) + XX;
    end
end

function [H, Qtsa, Msi, Tsi, R] = extractsol(x)

    global n m nun aux

    H    = zeros(n,m);
    Qtsa = zeros(n,m);
    Msi  = zeros(n,m);
    Tsi  = zeros(n,m);

    dim = n*m*nun;
    len = n*m*nun+aux;

    H(:)    = x(1:nun:dim);
    Qtsa(:) = x(2:nun:dim);
    Msi(:)  = x(3:nun:dim);
    Tsi(:)  = x(4:nun:dim);

    if (aux > 0)
        R = zeros(aux,1);
        R = x(dim+1:len);
    else
        R = [];
    end
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