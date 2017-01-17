function wrtbcsr (beg,jco,co, fnm, dir)
% WRTBCSR  writes a matrix in CSR format to a FORTRAN unformatted file.
% WRTBCSR(BEG,JCO,CO,FNM)     Writes the matrix stored in CSR format
%                             to the FORTRAN unformatted file FNM in the
%                             current directory.
% WRTBCSR(BEG,JCO,CO,FNM,DIR) Writes the matrix stored in CSR format
%                             to the FORTRAN unformatted file ~/DIR/FNM.
% 
%  wrtbcsr (beg,jco,co,'SA','/matlab/test/')
%  writes the matrix in CSR format (beg, jco, co) to the file
%  ~/matlab/test/SA
%
%  The file to be created will consist of 4 consecutive blocks.
%  Each block consists of
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
if nargin ~= 4 
 % dir=[getenv('HOME'),'/',dir];
  if dir(length(dir)) ~= '/'
    dir = [dir, '/'];
  end
  fnm = [dir,fnm];
end;
%
fid = fopen (fnm, 'wb');
if fid == -1
  error (sprintf('Cannot open file:  %s', fnm))
end
%
n   = length(beg) - 1;
nnz = beg(n+1) - 1;
%
bc  = 4*1;
cnt = fwrite(fid, bc, 'int32');   
cnt = fwrite(fid, n , 'int32');
cnt = fwrite(fid, bc, 'int32');
%
bc  = 4*(n+1);
cnt = fwrite(fid, bc , 'int32');
cnt = fwrite(fid, beg, 'int32');
cnt = fwrite(fid, bc , 'int32');
%
bc  = 4*nnz;
cnt = fwrite(fid, bc , 'int32');
cnt = fwrite(fid, jco, 'int32');
cnt = fwrite(fid, bc , 'int32');
%
bc  = 8*nnz;
cnt = fwrite(fid, bc , 'int32');
cnt = fwrite(fid, co , 'float64');
cnt = fwrite(fid, bc , 'int32');
%  
fclose(fid);
%disp(sprintf('CPU time (wrtbcsr): %g', cputime-tt))
