m = 16;
N = m*m;
ksub = 16;
ksup = 16;
ldimAB = 2*ksub+ksup+1;
AB = zeros(ldimAB, N);
kdiag = ksub + 1 + ksup;

for icol = 1:N
    i1 = max(1, icol-ksup);
    i2 = min(N, icol+ksub);
    for irow = i1:i2
        irowb = irow - icol + kdiag;
        AB(irowb, icol) = A(irow, icol);
    end
end
    