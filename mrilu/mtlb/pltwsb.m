function pltwsb(fnm)
% PLTWSB  Plot the allocation of segments in the WS buffer.
% 
% PLTWSB(FNM)  Plots the layout of the claimed and free segments in
%              the Work Space buffer.  The claimed segments are indicated
%              by a red area in the plot.
%              The file FNM  contains the index-value pairs where the
%              segments start and end, in the order of increasing indices.
%              The value is one for a claimed block and zero for a free
%              block.
%              The file FNM can be produced with the help of the Fortran
%              subroutine  wrtwsb.
%
if nargin ~= 1
  error ('Use a one input argument')
end
%
fid = fopen (fnm, 'rt');
if fid == -1
  error (sprintf('Cannot open file:  %s', fnm))
end
%
[bv,nr] = fscanf (fid, '%d');
% 
fclose(fid);
%
if nr < 4  |  floor(nr/4)*4 ~= nr
  error (sprintf('Number of integers in file, %i, is not correct!', nr))
end
%
bv    = reshape (bv, 2, floor(nr/2));
n     = nr/2;
polyg = [ bv, [bv(1,n); 0], [ bv(1,1); 0] ];
%
fill(polyg(1,:), polyg(2,:), 'r');
hold on
plot(bv(1,:), bv(2,:), 'y')
axis([bv(1,1) bv(1,n) -0.5 1.5]);
set (gca, 'YTick', [], 'YTickLabel', [])
title (sprintf('Allocation WS buffer, from %s ', fnm))
%
zoom on;
hold off
%
