C = load('atmos_jac'); C = spconvert(C);

B = load('atmos_B');

n = 6; m = 6; l = 1; dof = 2; 

dim = n*m*l*dof;

B = spdiags(B', 0, dim, dim);

idx = [];

for i = 1:dof;
    idx = [idx, i:dof:dim];
end

C   = C(idx,idx); % reordering

figure(1); 
spy(C);

figure(2);
[V,D]=eigs(C,B,10,0,opts);
plot(real(diag(D)),imag(diag(D)),'*')