cdata = dlmread('cdata.txt', ' ', 1, 0)

plot(cdata(:,1),cdata(:,2)); hold on
plot(cdata(:,1),cdata(:,3)); 
plot(cdata(:,1),cdata(:,4)); 
plot(cdata(:,1),cdata(:,5)); 
plot(cdata(:,1),cdata(:,6)); 

hold off

size(cdata)

xlabel('par')
grid on