[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');

otemp  = importdata('atmos_oceanTemp.txt');
state  = importdata('atmos_state.txt');

T0  = 15.0       ; %//! reference temperature
RtD = 180/pi;
To  = reshape(T0 + otemp,n,m);
Ta  = reshape(T0 + state,n,m);

figure(5)
contourf(RtD*x, RtD*y,To',15);
colorbar
title('SST ')
xlabel('Longitude')
ylabel('Latitude')
exportfig('oceanTemp.eps')

figure(6)
contourf(RtD*x,RtD*y,Ta',15);
colorbar
title('Atmosphere')
xlabel('Longitude')
ylabel('Latitude')
exportfig('atmosTemp.eps')
