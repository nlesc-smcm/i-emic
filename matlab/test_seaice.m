!killall vsm
J   = load('seaIceJac'); J = spconvert(J);
Jn  = load_numjac('seaIceNumJac');


n       = 4;
m       = 4;
dfa     = 4;
aux     = 1;
nseaice = n*m*dfa;
len     = nseaice + aux;

si_idx = [];
for i = 1:dfa
    si_idx = [si_idx, (i:dfa:nseaice)];
end
if (aux > 0)
    si_idx = [si_idx, si_idx(end)+1:len];
end


Jo  = J(si_idx,si_idx);
Jno = Jn(si_idx,si_idx);
