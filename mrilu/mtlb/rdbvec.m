function vec = rdbvec (fnm, dir)
% RDBVEC  reads a vector from unformatted file.
% VEC = RDBVEC(FNM)     Returns the vector stored in the unformatted file
%                       FNM in the current directory in the variabel VEC.
% VEC = RDBVEC(FNM,DIR) Returns the vector stored in the unformatted file
%                       FNM in the directory DIR (relative pathname with
%                       respect to the home directory) in the variable VEC. 
% 
%  vec = rdbvec ('RHS','/matlab/test/')
%  creates a vector variable from the file ~/matlab/test/RHS
%
%  The file has been created by a FORTRAN program or by the matlab function
%  WRTBVEC and consists of 2 consecutive blocks.  Each block consists of
%    4 byte integer     byte count of the data block
%    Actual data block
%    4 byte integer     byte count of the data block
%  
%  Data block 1:   1 integer N :  the length of a vector x
%  Data block 2:   N doubles:     the array  x  with values of the elements
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
if (bc1 ~= 8*n)
  error(inconsdata)
end
vec = fread(fid, n, 'float64');
bc2 = fread(fid, 1  , 'int32');
if (bc2 ~= bc1)
  error(inconsdata)
end
%  
fclose(fid);
%disp(sprintf('CPU time (rdbvec): %g', cputime-tt))
