!killall vsm&
clear all
JnC = load_numjac('JnC');

vsm(JnC)

[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = ...
    readfort44('fort.44');

dfo = 6; dfa = 2; aux = 1; dfs = 4;

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

%oce_idx = sort(oce_idx);

for i = 1:dfa
    atm_idx = [atm_idx, (nocean+i:dfa:nocean+natmos)];
end


for i = 1:aux
    atm_idx = [atm_idx, atm_idx(end) + i];
end

%atm_idx = sort(atm_idx);

for i = 1:dfs
    sei_idx = [sei_idx, (nocean+natmos+aux+i:dfs:nocean+natmos+aux+nseaice)];
end
%sei_idx = sort(sei_idx);

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
C12 = padding(C12, nocean, natmos + aux);
C13 = load('C_Ocean-SeaIce');       C13 = spconvert(C13);
C13 = padding(C13, nocean, nseaice);

C21 = load('C_Atmosphere-Ocean');   C21 = spconvert(C21);
C21 = padding(C21, natmos + aux, nocean);
C22 = load('J_Atmosphere');         C22 = spconvert(C22);
C23 = load('C_Atmosphere-SeaIce');  C23 = spconvert(C23);
C23 = padding(C23, natmos + aux, nseaice);

C31 = load('C_SeaIce-Ocean');       C31 = spconvert(C31);
C31 = padding(C31, nseaice, nocean);
    
C32 = load('C_SeaIce-Atmosphere');  C32 = spconvert(C32);
C32 = padding(C32, nseaice, natmos + aux);

C33 = load('J_SeaIce');             C33 = spconvert(C33);

% block might be smaller
tr12 = size(C12,1);
tr32 = size(C32,2);

C11 = C11(oce_idx, oce_idx);
C12 = C12(oce_idx, atm_idx); 
C13 = C13(oce_idx, sei_idx);

C21 = C21(atm_idx, oce_idx);
C22 = C22(atm_idx, atm_idx);
C23 = C23(atm_idx, sei_idx);

C31 = C31(sei_idx, oce_idx);
C32 = C32(sei_idx, atm_idx);
C33 = C33(sei_idx, sei_idx);

C = [C11, C12, C13; C21, C22, C23; C31, C32, C33];

