!killall vsm
!../../build/src/tests/test_atmos
old = load('Jno.mat')
J   = load('atmosJac'); J = spconvert(J);
Jn  = load_numjac('atmosNumJac');


n       = 8;
m       = 8;
dfa     = 3;
aux     = 1;
natmos  = n*m*dfa;

atm_idx = [];
for i = 1:dfa
    atm_idx = [atm_idx, (i:dfa:natmos)];
end

for i = 1:aux
    atm_idx = [atm_idx, atm_idx(end) + i];
end

Jo  = J(atm_idx,atm_idx);
Jno = Jn(atm_idx,atm_idx);
