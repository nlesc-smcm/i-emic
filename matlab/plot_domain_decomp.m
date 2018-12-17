function [] = plot_domain_decomp(n,m,nprocs)

    img = zeros(n,m);

    % using Jonas' dom decomp
    t1  = nprocs;
    t2  = 1;
    
    npM = t1;
    npN = t2;

    r = 0; % remainder
    r_min = 100;

    while (t1 > 0)
        t2 = floor(nprocs/t1);
        r  = abs(floor(m/t1) - floor(n/t2));

        if ((t1*t2 == nprocs) && (r <= r_min))
            r_min = r;
            npM   = t1;
            npN   = t2;
        end
        t1 = t1 - 1;
    end

    stepN = floor(n/npN);
    stepM = floor(m/npM);
    
    img(1:stepN:end,:)     = 1;
    img(stepN:stepN:end,:) = 1;
    img(:,1:stepM:end)     = 1;
    img(:,stepM:stepM:end) = 1;

    contour(img',1); 
    
end