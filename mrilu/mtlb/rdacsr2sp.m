 function A=rdacsr2sp(fnm, dir);
 if (nargin==1) dir=pwd; end
 [beg,jco,co] = rdacsr (fnm, dir);
 A = csr2sp (beg,jco,co);