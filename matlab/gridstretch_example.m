%--------------------------------------------
%  Grid stretching example
%  Ideally we have half of the grid points in
%  the upper 500m.
%--------------------------------------------

zmax  =     0;
zmin  =    -1;
L     =    12;
H     =  5000;

dz   = (zmax - zmin) / L;

ze  = zeros(L,1);
zwe = zeros(L,1);

for k = 1:L
  ze(k)  = (k-.5)*dz + zmin;
  zwe(k) = k*dz + zmin;
end

plot(ze*H,'.-'); hold on

qz = 1.8;
tz = fz(ze,qz); % transformed grid

plot(tz*H,'.-'); hold off

tz*H
grid on
