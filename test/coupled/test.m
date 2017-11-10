JnC = load_numjac('JnC');

n = 6; m = 6; l = 4; dfo = 6; dfa = 2;

surfb = find_row(dfo, n, m ,l, 1, 1, 4, 1);
surfe = find_row(dfo, n, m ,l, 6, 6, 4, dfo);

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

JnC22 = JnC(atm_idx, atm_idx);
JnC11 = JnC(oce_idx, oce_idx);
JnC12 = JnC(oce_idx, atm_idx);
JnC21 = JnC(atm_idx, oce_idx);

figure(1)
spy([JnC11,JnC12;JnC21,JnC22]);

figure(2)
spy(JnC)

C11 = load('C11'); C11 = spconvert(C11);
C12 = load('C12'); C12 = spconvert(C12);
C21 = load('C21'); C21 = spconvert(C21);
C22 = load('C22'); C22 = spconvert(C22);

aatm_idx = atm_idx - nocean;

figure(3)
spy(C12(surfb:surfe,:))



numC12 = JnC(surfb:surfe, surfe+1:end);
numC21 = JnC(surfe+1:end, surfb:surfe);

