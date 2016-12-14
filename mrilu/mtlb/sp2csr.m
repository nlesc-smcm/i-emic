function [beg,jco,co] = sp2csr (sa)
% SP2CSR   converts a MATLAB sparse matrix into CSR format.
% [BEG,JCO,CO] = SP2CSR(SA)  Returns the sparse matrix SA in the
%                            parameters BEG, JCO and CO (CSR format).
%  

%tt = cputime;
[m,n] = size(sa);
[jco,ico,co] = find(sa');
beg = cumsum(full([1,sum(sparse(jco,ico,1,m,n))]));
%
%disp(sprintf('CPU time (sp2csr): %g', cputime-tt))
