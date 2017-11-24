C = load('atmos_jac'); C = spconvert(C);

JnC = load_numjac('JnC');

n = 8; m = 8; l = 1; dof = 2; 

dim = n*m*l*dof;

figure(1); 
spy(C);
figure(2); 
spy(JnC);
figure(3); 
spy(abs(JnC-C)./abs(C)>1e-4);