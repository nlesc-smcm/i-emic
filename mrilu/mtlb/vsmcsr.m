function vsmcsr (beg,jco,co,fnm,dir)
% VSMCSR  Visualize Sparse Matrix given in CSR format.
% VSMCSR(BEG,JCO,CO)         Uses a temporary binary file to be read
%                            by  vsm  which visualizes the matrix given in
%                            CSR format (BEG,JCO,CO).
% VSMCSR(BEG,JCO,CO,FNM)     Creates the binary file FNM in the current
%                            directory to be used by  vsm.
% VSMCSR(BEG,JCO,CO,FNM,DIR) Creates the binary file ~/DIR/FNM to be used
%                            by  vsm.
%
vsmopt = [' '];
if nargin==3
  vsmopt = ['-d', vsmopt];
  fnm = tempname;
elseif nargin ~= 4
  dir=[getenv('HOME'),'/',dir];
  if dir(length(dir)) ~= '/'
    dir = [dir, '/'];
  end
  fnm = [dir,fnm];
end;
%
wrtbcsr(beg,jco,co,fnm);
%
eval(['!vsm ', vsmopt, fnm , ' &']);
