function [] = exportfig(string,fs,dims)
if nargin < 3
    dims = [14,10];
end
if nargin < 2
    fs = 10;
end

width = dims(1); height = dims(2);
set(gca,'fontsize',fs)
set(gcf,'paperunits','centimeters')
set(gcf,'paperposition',[1,1,width,height])
saveas(gcf, string, 'epsc2')
