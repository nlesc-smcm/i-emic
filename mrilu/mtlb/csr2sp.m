function SA = csr2sp (beg,jco,co)
% CSR2SP   converts a matrix in CSR format in a MATLAB sparse array.
% CSR2SP(BEG,JCO,CO)  Returns the sparse matrix stored in CSR format in the
%                     parameters BEG, JCO and CO.
%  

%tt = cputime;
n   = length(beg) - 1;
nnz = beg(n+1) - 1;
ico = zeros(nnz,1);
for i=1:n
  ico(beg(i):beg(i+1)-1) = i*ones(beg(i+1)-beg(i),1);
end;
%
% SA = sparse(ico,jco,co);  % More time consuming then next 2 statements!!
%  
SA = sparse(jco,ico,co);
SA = SA';
%disp(sprintf('CPU time (csr2sp): %g', cputime-tt))
