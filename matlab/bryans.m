wfun = @(yy) 0.2 - 0.8*sin(6*abs(yy)) - 0.5*(1 - tanh(10*abs(yy)))...
    -0.5*(1 - tanh(10*(pi/2 - abs(yy))));

ymin = 10; %degrees
ymax = 80;  
res  = 0.1;
yy = (ymin:res:ymax)*pi/180;

plot(wfun(yy), yy*180/pi,'linewidth',2)

ylabel('latitude')
xlabel('wind forcing')

exportfig('bryans.eps',10,[5,10]);