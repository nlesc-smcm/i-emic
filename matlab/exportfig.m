function [] = exportfig(string,fs,dims,invert)
    if nargin < 4
        invert = false;
    end
    if nargin < 3
        dims = [14,10];
    end
    if nargin < 2
        fs = 10;
    end

    width = dims(1); height = dims(2);
    set(gca,'fontsize',fs)
    set(gcf,'color','w')

    set(gcf,'paperunits','centimeters')
    set(gcf,'paperposition',[1,1,width,height])
    set(gcf,'papersize',[width,height])

    if invert
        invertcolors();        
    end
    
    set(gcf,'inverthardcopy','off')
    fprintf('saving to %s...\n', string);
    %print(string, '-depsc2', '-painters', '-opengl');
    print(string, '-depsc2', '-painters');
    %print(string, '-depsc2');
    %print(string)
    %saveas(gcf,string);
end

