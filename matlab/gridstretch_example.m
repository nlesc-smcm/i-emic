%--------------------------------------------
%  Grid stretching example
%--------------------------------------------

zmax = 0;
zmin = -1;
L = 16;

dz = (zmax - zmin) / L;

ze  = zeros(L,1);
zwe = zeros(L,1);

for k = 1:L
	ze(k)  = (k-.5)*dz + zmin
	zwe(k) = k*dz + zmin;
end

plot(ze); hold on

qz = 2.3;

plot(fz(ze,qz)); hold off
