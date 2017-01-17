function wrtavec (x, fnm, dir)
% WRTAVEC    Writes a vector to ASCII file.
% WRTAVEC(X,FNM)     Writes the vector X, length and values, to the
%                    ASCII file FNM in the current directory.
% WRTAVEC(X,FNM,DIR) Writes the vector X, length and values, to the
%                    ASCII file ~/DIR/FNM.
% 
%  wrtavec (x,'SA','/matlab/test/')
%  writes the vector in format (n, x) to the file
%  ~/matlab/test/SA
%
%
% One number per line is written to the file.
% The order of the numbers in the file are
%    integer:              the length of the vector N
%    N * double precision: the values of the vector elements.
%

%tt = cputime;
if nargin ~= 2 
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
n   = length(x);
%
cnt = fprintf(fid, '%d\n', n);
%
cnt = fprintf(fid, '%+22.15e\n', x);
%  
fclose(fid);
%disp(sprintf('CPU time (wrtavec): %g', cputime-tt))
