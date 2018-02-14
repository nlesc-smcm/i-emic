    c1 = 3.8e-3;  
    c2 = 21.87;   
    c3 = 265.5;   
    c4 = 17.67;   
    c5 = 243.5;   

    t0  = 15;
    ti0 = 0;
    q0  = 8e-3;
    dq  = 1e-3; % typical deviation
    dt  = 1;    % typical deviation

    rhoa_    = 1.25;
    rhoo_    = 1024;
    hdima_   = 8400.;
    hdimq_   = 1800.;
    hdim_    = 4000.;
    cpa_     = 1000.;
    d0_      = 3.1e+06;
    kappa_   = 1e+06;
    arad_    = 212.0;
    brad_    = 1.5;
    sun0_    = 1360.;
    c0_      = 0.43;
    ce_      = 1.3e-03;
    ch_      = 0.94 * ce_;
    uw_      = 8.5;
    t0a_     = 15.0; 
    t0o_     = 15.0;  
    t0i_     = 0.0;   
    tdim_    = 1.0; 
    q0_      = 15e-3; 
    qdim_    = 0.01;  
    lv_      = 2.5e06; 
    
    udim_    = 0.1;
    r0dim_   = 6.37e06;

    eta_  = (rhoa_ / rhoo_) * ce_ * uw_;
    muoa_ =  rhoa_ * ch_ * cpa_ * uw_;
    
    % background saturation humidity
    qso  = @(T) c1*exp(c4*T/(T+c5));
    qsi  = @(T) c1*exp(c2*T/(T+c3));
    qso0 = qso(t0)
    qsi0 = qsi(ti0)

    w = 1e-6;

    E = @(T,q) rhoa_ * ce_ * uw_ / rhoo_ * ( qso(T) - q);

    Etyp = E(t0, q0) + (E(t0+w,q0)-E(t0,q0))/w * dt + (E(t0,q0+w)-E(t0,q0))/w ...
           * dq

    E0 = E(t0,q0);
    P0 = E0
    
    % E(t0+dt, q0+dq)
    % 
    % dqso0 = (qso(t0+1e-6)-qso(t0))/1e-6
    % 
    % dqso0 = (c1 * c4 * c5) / (t0 + c5)^2 *exp( (c4 * t0) / (t0 + c5) ...
    %                                            ) 
    % eta_  = (rhoa_ / rhoo_) * ce_ * uw_;
    %
    %    (E(t0+w,q0)-E(t0,q0))/w                                        
    %    
    %    dqso0 * eta_
    %    
    %    rhoo_ * lv_ * Etyp
    %    
    %    (E(t0,q0+w)-E(t0,q0))/w
    %    
    %    -eta_
        
    