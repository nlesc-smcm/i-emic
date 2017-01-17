function [vec] = rdavec (fnm, dir)
% RDAVEC  reads a vector from an ASCII file.
% [VEC] = RDAVEC(FNM)     Returns the vector stored in the ASCII file
%                         FNM in the current directory.
% [VEC] = RDAVEC(FNM,DIR) Returns the vector stored in the ASCII file
%                         FNM in the directory  DIR (relative pathname
%                         with respect to the home directory). 
% 
%  [vec] = rdavec ('SA','/matlab/test/')
%  creates a vector from the file ~/matlab/test/SA
%
%  The file consists of 2 consecutive blocks.  Each block starts on a
%  new line.
%  Data block 1:   1 integer N :  the length of the vector.
%  Data block 2:   N doubles:     the values of the elements of the
%                                 vector.
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
fid = fopen (fnm, 'r');
if fid == -1
  error (sprintf('Cannot open file:  %s', fnm))
end
%
n   = fscanf(fid, '%d', 1);
%
vec = fscanf(fid, '%e', n);
%
% 
fclose(fid);
%disp(sprintf('CPU time (rdavec): %g', cputime-tt))
