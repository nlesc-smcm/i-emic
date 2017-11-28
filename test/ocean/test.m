C = load('ocean_jac'); C = spconvert(C);

JnC = load_numjac('ocean_numjac');

n = 6; m = 6; l = 4; dof = 6; 

dim = n*m*l*dof;

idx = [];

for i = 1:dof;
    idx = [idx, i:dof:dim];
end

C   = C(idx,idx); % reordering
JnC = JnC(idx,idx); % reordering

figure(1); 
spy(C);

figure(2); 
spy(JnC);

figure(3); 

spy(abs(JnC-C)./abs(C)>1e-4);