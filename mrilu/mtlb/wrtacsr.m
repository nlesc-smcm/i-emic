function wrtacsr (beg,jco,co, fnm, dir)
% WRTACSR     Writes sparse matrix, in CSR format, to ASCII file.
% WRTACSR(BEG,JCO,CO,FNM)     Writes the matrix stored in CSR format
%                             to the FORTRAN unformatted file FNM in the
%                             current directory.
% WRTACSR(BEG,JCO,CO,FNM,DIR) Writes the matrix stored in CSR format
%                             to the FORTRAN unformatted file ~/DIR/FNM.
% 
%  wrtacsr (beg,jco,co,'SA','/matlab/test/')
%  writes the matrix in CSR format (beg, jco, co) to the file
%  ~/matlab/test/SA
%
% One number per line is written to the file.
% The order of the numbers in the file are
%    integer:        the order of the matrix N
%    N+1 * integer:  indices indicating the beginning of a row
%    NNZ * integer:  column numbers, with NNZ := begA(N+1)-1
%    NNZ * double precision: nonzero values
%

%tt = cputime;
if nargin ~= 4 
 % dir=[getenv('HOME'),'/',dir];
  if dir(length(dir)) ~= '/'
    dir = [dir, '/'];
  end
  fnm = [dir,fnm];
end;
%
fid = fopen (fnm, 'wt');
if fid == -1
  error (sprintf('Cannot open file:  %s', fnm))
end
%
n   = length(beg) - 1;
nnz = beg(n+1) - 1;
%
cnt = fprintf(fid, '%d\n', n);
%
cnt = fprintf(fid, '%d\n', beg);
%
cnt = fprintf(fid, '%d\n', jco);
%
cnt = fprintf(fid, '%+22.15e\n', co);
%  
fclose(fid);
%disp(sprintf('CPU time (wrtacsr): %g', cputime-tt))
