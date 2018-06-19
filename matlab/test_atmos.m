J  = load('atmosJac'); J = spconvert(J);
Jn = load_numjac('atmosNumJac');


n   = 8
m   = 8
dfa = 3
natmos  = n*m*dfa;

assert(size(Jn,1) == natmos)

atm_idx = [];
for i = 1:dfa
    atm_idx = [atm_idx, (i:dfa:natmos)];
end

vsm(Jn(atm_idx,atm_idx))
