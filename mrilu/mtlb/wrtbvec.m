function wrtbvec (x, fnm, dir)
% WRTBVEC  writes a vector to a FORTRAN unformatted file.
% WRTBVEC(X,FNM)     Writes the vector X, length and values, to the
%                    FORTRAN unformatted file FNM in the current directory.
% WRTBVEC(X,FNM,DIR) Writes the vector X, length and values, to the
%                    FORTRAN unformatted file ~/DIR/FNM.
% 
%  wrtbvec (x,'SA','/matlab/test/')
%  writes the vector in format (n, x) to the file
%  ~/matlab/test/SA
%
%  The file to be created will consist of 2 consecutive blocks.
%  Each block consists of
%    4 byte integer     byte count of the data block
%    Actual data block
%    4 byte integer     byte count of the data block
%  
%  Data block 1:   1 integer N :  the length of the vector
%  Data block 2:   N doubles   :  the array  x
%

%tt = cputime;
if nargin ~= 2 
  dir=[getenv('HOME'),'/',dir];
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
n   = length(x);
%
bc  = 4*1;
cnt = fwrite(fid, bc, 'int32');   
cnt = fwrite(fid, n , 'int32');
cnt = fwrite(fid, bc, 'int32');
%
bc  = 8*n;
cnt = fwrite(fid, bc , 'int32');
cnt = fwrite(fid, x,   'float64');
cnt = fwrite(fid, bc , 'int32');
%  
fclose(fid);
%disp(sprintf('CPU time (wrtbvec): %g', cputime-tt))
