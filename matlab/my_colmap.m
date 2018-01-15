function [out] = my_colmap(colrange, mn, sze)

% colrange : color range (caxis)
% mn       : data mean
% sze      : number of colors

    if nargin < 3
        sze = 128;
    end

    if nargin < 2
        mn = mean(colrange);
    end

    cmax = max(colrange);
    cmin = min(colrange);
    clen = abs(cmax - cmin);

    if cmin > mn
        N1 = 0;
    else
        N1 = round(abs(cmin-mn)/clen * sze);
    end

    if cmax < mn
        N2 = 0;
    else
        N2 = round(abs(cmax-mn)/clen * sze);
    end
    
    par = [0    0.4470    0.7410;  0.8500    0.3250    0.0980];
    neg  = par(1,:);
    pos  = par(2,:);
    mid  = [1,1,1];
    col1 = [linspace(neg(1),mid(1),N1)',linspace(neg(2),mid(2),N1)',linspace(neg(3),mid(3),N1)'];
    col2 = [linspace(mid(1),pos(1),N2)',linspace(mid(2),pos(2),N2)',linspace(mid(3),pos(3),N2)'];
    out  = [col1;col2];
end