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

figure(11); 
spy(C);

figure(12); 
spy(JnC);

figure(13); 
diff11 = abs(C - JnC)./abs(C) > tol;
if sum(diff11(:)) > 0
    figure(13)
    diff11 = (abs(C - JnC));
    imagesc(diff11); colorbar
    title('diff')
end

