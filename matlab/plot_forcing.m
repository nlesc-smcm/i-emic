maskfile='fort.44';

[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = ...
readfort44(maskfile);

%% - DEFINE CONSTANTS - ----------------------------------------------
udim  = 0.1;                 %[m/s]   Velocity scale
r0dim = 6.4e6;               %[m]     Radius of Earth
T0    = 15;                  %[deg C] Reference temperature
S0    = 35;                  %[psu]   Reference salinity
RtD   = 180/pi;              %[-]     Radians to degrees

%% - F.Bryans (1987) idealized wind profile  
wfun = @(yy) 0.2 - 0.8*sin(6*abs(yy)) - 0.5*(1 - tanh(10*abs(yy)))...
       -0.5*(1 - tanh(10*(pi/2 - abs(yy))));

%% - Idealized temperature/salinity profile
temfun = @(yy) cos(pi*(yy) / (ymax));

figure(1)
Y = linspace(ymin,ymax,100);
plot(Y*RtD, wfun(Y),'color',[0,0,.85],'linewidth',2)
xlabel('Latitude')
grid on
hold on
plot(Y*RtD, temfun(Y), 'color',[.85,0,0], 'linewidth',2)
hold off
ylim([-1.1,1.1]);
xlim([-80,80]);
%legend('\tau^\phi','T_S, F_S');
exportfig('forcing.eps',11,[15,6]);


