function [] = invertcolors()
    set(gca,'color','k');
    set(gca,'xcolor','w');
    set(gca,'ycolor','w');
    set(gca,'gridcolormode','manual');
    set(gca,'minorgridcolormode','manual');
    set(gca,'gridcolor',[.15,.15,.15]);
    set(gca,'minorgridcolor',[.15,.15,.15]);
    set(gca,'gridalpha',1);
    set(gca,'minorgridalpha',1);

    set(gcf,'color','k');
    set(legend,'color','k');
    set(legend,'edgecolor',[.5,.5,.5]);
    set(legend,'textcolor','w');        

    title = get(gca,'title');
    title.Color=[1,1,1];
end