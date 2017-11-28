clear all
JnC = load_numjac('JnC');

n = 6; m = 6; l = 4; dfo = 6; dfa = 2; aux = 1

surfb = find_row(dfo, n, m ,l, 1, 1, l, 1);
surfe = find_row(dfo, n, m ,l, n, m, l, dfo);

nocean = n*m*l*dfo;
natmos = n*m*dfa;

atm_idx = [];
oce_idx = [];

for i = 1:dfo
    oce_idx = [oce_idx, (i:dfo:nocean)];
end

for i = 1:dfa
    atm_idx = [atm_idx, (nocean+i:dfa:nocean+natmos)];
end

for i = 1:aux
    atm_idx = [atm_idx, atm_idx(end) + i];
end

JnC22 = JnC(atm_idx, atm_idx);
JnC11 = JnC(oce_idx, oce_idx);
JnC12 = JnC(oce_idx, atm_idx);
JnC21 = JnC(atm_idx, oce_idx);

C11 = load('C11'); C11 = spconvert(C11);
C12 = load('C12'); C12 = spconvert(C12);
C21 = load('C21'); C21 = spconvert(C21);
C22 = load('C22'); C22 = spconvert(C22);

numC12 = JnC(surfb:surfe, surfe+1:end);
numC21 = JnC(surfe+1:end, surfb:surfe);
numC11 = JnC(1:surfe, 1:surfe);
numC22 = JnC(surfe+1:end, surfe+1:end);

numC11 = numC11(oce_idx,oce_idx);
C11    =    C11(oce_idx,oce_idx);

atm_idx = atm_idx - surfe;
numC22 = numC22(atm_idx,atm_idx);
C22    =    C22(atm_idx,atm_idx);

figure(1)
spy(C11)
title('C11')


figure(2)
spy(JnC11)
title('num C11')

tol = 1e-4

diff11 = abs(C11 - numC11)./abs(C11) > tol;
if sum(diff11(:)) > 0
    figure(3) 
    spy(diff11)
    figure(4)
    diff11 = abs(C11 - numC11)
    imagesc(diff11); colorbar
    axis square
    title('C11 diff')
end

return 

C12    = C12(surfb:end,:);
numC12 = numC12(1:end-1,:);
diff12 = abs(C12 - numC12)./abs(C12) > tol;
if sum(diff12(:)) > 0
    figure(4)
    spy(diff12)
    title('C12 diff')
end

C21    = C21(:,surfb:end);
numC21 = numC21(:,1:end-1);
diff21 = abs(C21 - numC21)./abs(C21) > tol;
if sum(diff21(:)) > 0
    figure(5)
    spy(diff21)
    title('C21 diff')    
end

diff22 = abs(C22 - numC22)./abs(C22) > tol;
if sum(diff22(:)) > 0
    figure(6)
    spy(diff22)
    title('C22 diff')    
end


