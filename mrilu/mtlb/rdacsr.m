function [beg,jco,co] = rdacsr (fnm, dir)
% RDACSR  reads a matrix in CSR format from an ASCII file.
% [BEG,JCO,CO] = RDACSR(FNM)     Returns the matrix stored in CSR format
%                                in the ASCII file FNM in the current
%                                directory.
% [BEG,JCO,CO] = RDACSR(FNM,DIR) Returns the matrix stored in CSR format
%                                in the ASCII file FNM in the directory
%                                DIR (relative pathname with respect to
%                                the home directory). 
% 
%  [beg,jco,co] = rdacsr ('SA','/matlab/test/')
%  creates a matrix in CSR format from the file ~/matlab/test/SA
%
%  The file consists of 4 consecutive blocks.  Each block starts on a
%  new line.
%  Data block 1:   1 integer N :  the row dimension of the matrix
%  Data block 2:   N+1 integers:  the array  beg  with pointers in jco and co
%                                 NNZ := beg(N+1)-1
%  Data block 3:   NNZ integers:  the array  jco  with row indices
%  Data block 4:   NNZ doubles:   the array  co  with values of the elements
%

%tt = cputime;
if nargin ~= 1
% dir=[getenv('HOME'),'/',dir];
  if dir(length(dir)) ~= '/'
    dir = [dir, '/'];
  end
  fnm = [dir,fnm];
end;
%
fid = fopen (fnm, 'r');
if fid == -1
  error (sprintf('Cannot open file:  %s', fnm))
end
%
n   = fscanf(fid, '%d', 1);
%
beg = fscanf(fid, '%d', n+1);
%
nnz = beg(n+1)-1;
jco = fscanf(fid, '%d', nnz);
%
co  = fscanf(fid, '%e', nnz);
% 
fclose(fid);
%disp(sprintf('CPU time (rdacsr): %g', cputime-tt))
