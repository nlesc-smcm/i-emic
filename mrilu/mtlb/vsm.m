function vsm(A,fnm,dir)
% VSM  Visualize Sparse Matrix
% VSM(A)         Uses a temporary binary file to be read by  vsm
%                which visualizes the sparse matrix A.
% VSM(A,FNM)     Creates the binary file FNM in the current directory
%                to be used by  vsm.
% VSM(A,FNM,DIR) Creates the binary file ~/DIR/FNM to be used
%                by  vsm.
%
vsmopt = [' '];
if nargin==1
    vsmopt = ['-d', vsmopt];
    fnm = tempname;
elseif nargin ~= 2
    dir=[getenv('HOME'),'/',dir];
    if dir(length(dir)) ~= '/'
        dir = [dir, '/'];
    end
    fnm = [dir,fnm];
end;
%
% Create the arrays beg, jco and co
% 
[jco,ico,co]=find(A');
beg=cumsum(full([1,sum(sparse(jco,ico,1))]));
%
wrtbcsr(beg,jco,co,fnm);
%
out = system(['vsm ', vsmopt, fnm , ' &']); 
fprintf(' vsm exit status: %d\n', out);

end