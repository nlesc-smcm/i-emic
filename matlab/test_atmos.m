C = load('atmos_jac'); C = spconvert(C);

B = load('atmos_B');

dim = size(C,1);
B = spdiags(B', 0, dim, dim);

figure(1); 
spy(C);

figure(2);
opts.maxit = 1000;
opts.tol = 1e-12;
neig = 6;
[V,D]=eigs(C,B,neig,0,opts);
for i = 1:neig
    fprintf('%3d: real: %2.5e imag: %2.5e \n', i, real(D(i,i)), ...
            imag(D(i,i)));
end
plot(real(diag(D)),imag(diag(D)),'*')