% Script voor het inlezen van de ijle matrix  S  en de factoren
% L, D  en  U  van de incomplete factorisatie van  S.
% Verder worden de volgende matrices berekend:
%   R   := S - L*D*U
%   RS  := R | NZ(S): de restrictie van R op de non-zeros van  S.
%   RLU := R | NZ(L+U): de restrictie van R op de non-zeros van  L+U.
%
[begS,jcoS,coS] = rdbcsr('LastSchurComp');
S = csr2sp(begS,jcoS,coS);
clear begS jcoS coS
[ms,ns]=size(S);
% 
[begL,jcoL,coL] = rdbcsr('ILDUfactor.L'); 
L = csr2sp(begL,jcoL,coL);
clear begL jcoL coL
[begU,jcoU,coU] = rdbcsr('ILDUfactor.U');
U = csr2sp(begU,jcoU,coU);
clear begU jcoU coU
[begD,jcoD,coD] = rdbcsr('ILDUfactor.D');
D = csr2sp(begD,jcoD,coD);
clear begD jcoD coD
%
R = S - L*D*U;
%
[is,js] = find(S);
nzs=length(is);
RS=zeros(nzs,1);
for k=1:nzs
   RS(k)=R(is(k),js(k));
end
RS = sparse(is,js,RS,ms,ns);
disp(sprintf('\nS is een %1ix%1i matrix met %1i non-zeros.', ms, ns, nzs))
clear is js nzs
% 
[ilu,jlu] = find(L+U);
nzlu=length(ilu);
RLU=zeros(nzlu,1);
for k=1:nzlu
   RLU(k)=R(ilu(k),jlu(k));
end
RLU = sparse(ilu,jlu,RLU,ms,ns);
%
disp(sprintf('Rest matrix R heeft %1i non-zeros.', nzlu))
clear ilu jlu nzlu
%
mrsR = max(abs(sum(R'))); maxrowsumR = mrsR(1,1);
disp(sprintf('Maximum absolute waarden rijsommen in R: %e\n', maxrowsumR))
clear mrsR
%
% End of  sldu.
