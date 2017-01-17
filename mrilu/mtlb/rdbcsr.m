function [beg,jco,co] = rdbcsr (fnm, dir)
% RDBCSR  reads a matrix in CSR format from unformatted file.
% [BEG,JCO,CO] = RDBCSR(FNM)     Returns the matrix stored in CSR format
%                                in the unformatted file FNM in the
%                                current directory.
% [BEG,JCO,CO] = RDBCSR(FNM,DIR) Returns the matrix stored in CSR format
%                                in the unformatted file FNM in the
%                                directory DIR (relative pathname with
%                                respect to the home directory). 
% 
%  [beg,jco,co] = rdbcsr ('SA','/matlab/test/')
%  creates a matrix in CSR format from the file ~/matlab/test/SA
%
%  The file has been created by a FORTRAN program or by the matlab function
%  WRTBCSR and consists of 4 consecutive blocks.  Each block consists of
%    4 byte integer     byte count of the data block
%    Actual data block
%    4 byte integer     byte count of the data block
%  
%  Data block 1:   1 integer N :  the row dimension of the matrix
%  Data block 2:   N+1 integers:  the array  beg  with pointers in jco and co
%                                 NNZ := beg(N+1)-1
%  Data block 3:   NNZ integers:  the array  jco  with row indices
%  Data block 4:   NNZ doubles:   the array  co  with values of the elements
%

%tt = cputime;
if nargin ~= 1 
  dir=[getenv('HOME'),'/',dir];
  if dir(length(dir)) ~= '/'
    dir = [dir, '/'];
  end
  fnm = [dir,fnm];
end;
%
fid = fopen (fnm, 'rb');
if fid == -1
  error (sprintf('Cannot open file:  %s', fnm))
end
inconsdata = sprintf('Inconsistent data in  %s', fnm);
%
bc1 = fread(fid, 1, 'int32');   
if (bc1 ~= 4)
  error(inconsdata)
end
n   = fread(fid, 1, 'int32');
bc2 = fread(fid, 1, 'int32');
if (bc2 ~= bc1)
  error(inconsdata)
end
%
bc1 = fread(fid, 1  , 'int32');
if (bc1 ~= 4*(n+1))
  error(inconsdata)
end
beg = fread(fid, n+1, 'int32');
bc2 = fread(fid, 1  , 'int32');
if (bc2 ~= bc1)
  error(inconsdata)
end
%
nnz = beg(n+1)-1;
bc1 = fread(fid, 1  , 'int32');
if (bc1 ~= 4*nnz)
  error(inconsdata)
end
jco = fread(fid, nnz, 'int32');
bc2 = fread(fid, 1  , 'int32');
if (bc2 ~= bc1)
  error(inconsdata)
end
%
bc1 = fread(fid, 1  , 'int32');
if (bc1 ~= 8*nnz)
  error(inconsdata)
end
co  = fread(fid, nnz, 'float64');
bc2 = fread(fid, 1  , 'int32');
if (bc2 ~= bc1)
  error(inconsdata)
end
%  
fclose(fid);
%disp(sprintf('CPU time (rdbcsr): %g', cputime-tt))
