function sp2csrwrta(sa,fnm,dir);
if (nargin == 2) dir=pwd; end
[beg,jco,co] = sp2csr (sa);
wrtacsr (beg,jco,co, fnm, dir);