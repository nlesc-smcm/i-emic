wfun = @(yy) 0.2 - 0.8*sin(6*abs(yy)) - 0.5*(1 - tanh(10*abs(yy)))...
    -0.5*(1 - tanh(10*(pi/2 - abs(yy))));

ymin = -80; %degrees
ymax = 80;  
res  = 0.1;
yy = (ymin:res:ymax)*pi/180;

plot(yy*180/pi, wfun(yy))