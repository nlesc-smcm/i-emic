[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');

otemp  = importdata('atmos_oceanTemp.txt');
state  = importdata('atmos_state.txt');

t0_     = 15.0       ; %//! reference temperature

To  = reshape(t0_ + otemp,n,m);
Ta  = reshape(t0_ + state,n,m);

figure(1)
contourf(x,y,To',15);
colorbar
title('SST')
xlabel('Longitude')
ylabel('Latitude')
exportfig('oceanTemp.eps')

figure(2)
contourf(x,y,Ta',15);
colorbar
title('Atmosphere')
xlabel('Longitude')
ylabel('Latitude')
exportfig('atmosTemp.eps')
