clear all
JnC = load_numjac('JnC');

n = 6; m = 6; l = 4; dfo = 6; dfa = 2; aux = 1; dfs = 3;

surfb = find_row(dfo, n, m ,l, 1, 1, l, 1);
surfe = find_row(dfo, n, m ,l, n, m, l, dfo);

nocean  = n*m*l*dfo;
natmos  = n*m*dfa;
nseaice = n*m*dfs;

ndim = nocean + natmos + aux + nseaice;

atm_idx = [];
oce_idx = [];
sei_idx = [];

for i = 1:dfo
    oce_idx = [oce_idx, (i:dfo:nocean)];
end

for i = 1:dfa
    atm_idx = [atm_idx, (nocean+i:dfa:nocean+natmos)];
end

for i = 1:aux
    atm_idx = [atm_idx, atm_idx(end) + i];
end

for i = 1:dfs
    sei_idx = [sei_idx, (nocean+natmos+aux+i:dfs:nocean+natmos+aux+nseaice)];
end

JnC11 = JnC(oce_idx, oce_idx);
JnC22 = JnC(atm_idx, atm_idx);
JnC33 = JnC(sei_idx, sei_idx);

JnC12 = JnC(oce_idx, atm_idx);
JnC13 = JnC(oce_idx, sei_idx);

JnC21 = JnC(atm_idx, oce_idx);
JnC23 = JnC(atm_idx, sei_idx);

JnC31 = JnC(sei_idx, oce_idx);
JnC32 = JnC(sei_idx, atm_idx);

tot_idx = [oce_idx, atm_idx, sei_idx];
JnC = JnC(tot_idx, tot_idx);

% remove offsets of ranges
oce_idx = oce_idx - 0;
atm_idx = atm_idx - nocean;
sei_idx = sei_idx - nocean - natmos - aux;

C11 = load('J_Ocean');              C11 = spconvert(C11);
C12 = load('C_Ocean-Atmosphere');   C12 = spconvert(C12);
%C13 = load('C_Ocean-SeaIce');       C13 = spconvert(C13);

C21 = load('C_Atmosphere-Ocean');   C21 = spconvert(C21);
C22 = load('J_Atmosphere');         C22 = spconvert(C22);
%C23 = load('C_Atmosphere-SeaIce');  C23 = spconvert(C23);

C31 = load('C_SeaIce-Ocean');       C31 = spconvert(C31);
C32 = load('C_SeaIce-Atmosphere');  C32 = spconvert(C32);
C33 = load('J_SeaIce');             C33 = spconvert(C33);

% block might be smaller
tr12 = size(C12,1);
tr32 = size(C32,2);

C11 = C11(oce_idx, oce_idx);
C12 = C12(oce_idx(1:tr12), atm_idx);
%C13 = C13(oce_idx, sei_idx);

C21 = C21(atm_idx, oce_idx(1:tr12));
C22 = C22(atm_idx, atm_idx);
%C23 = C23(atm_idx, sei_idx);

C31 = C31(sei_idx(1:size(C32,1)), oce_idx);
C32 = C32(sei_idx(1:size(C32,1)), atm_idx(1:size(C32,2)));
C33 = C33(sei_idx, sei_idx);

