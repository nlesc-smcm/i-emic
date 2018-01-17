C = load('ocean_jac'); C = spconvert(C);

%JnC = load_numjac('ocean_numjac');

B = load('ocean_B');

n = 8; m = 8; l = 8; dof = 6; 

dim = n*m*l*dof;

B = spdiags(B', 0, dim, dim);

%idx = [];
%
%for i = 1:dof;
    %    idx = [idx, i:dof:dim];
    %end

%C   = C(idx,idx); % reordering
%JnC = JnC(idx,idx); % reordering

figure(1); 
spy(C);

%figure(2); 
%spy(JnC);

figure(3); 
tic
opts.maxit = 1000;
opts.tol = 1e-12;
neig = 5;
[V,D]=eigs(C,B,neig,0,opts);
toc

for i = 1:neig
    fprintf('%3d: real: %2.5e imag: %2.5e \n', i, real(D(i,i)), ...
            imag(D(i,i)));
end

plot(real(diag(D)),imag(diag(D)),'*')