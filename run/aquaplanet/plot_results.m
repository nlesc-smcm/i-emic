basenames={'eps0.05_','eps0.1_','eps0.2_','eps0.4_'};
for i = 1:4
    basename = basenames{i};
    opts= [];
    cols   = lines(3);
    styles = {{'color',cols(1,:),'linewidth',1.5},...
              {'color',cols(2,:),'linewidth',1.5},...
              {'color',cols(3,:),'linewidth',1.5}};

    ctr = 0;

    for i = [8,16,32]
        ctr = ctr + 1;
        file=[basename,num2str(i),'/cdata.txt'];
        if i == 8
            opts.hold=false;
        else
            opts.hold=true;
        end

        opts.invert=false;
        opts.plot_entry=15;
        opts.point=-1;
        opts.lsty=styles{ctr};    
        [~,~,h(ctr)] = plot_cdata(file,opts);
        
    end
    xlim([0.5,1.6]);
    title('');
    ylabel('A^{si} / A');
    xlabel('Solar forcing \lambda_\Sigma');

    legend('m=8', 'm=16', 'm=32','location','southeast');

    h(1).Visible='off';
    h(2).Visible='off';
    h(3).Visible='off';  

    for i = 1:3
        h(i).Visible='on';
        fname = ['aquaplanbif',num2str(i),'.eps'];

        exportfig(fname, 13, [18,13], opts.invert);

        dname = ['/home/erik/Projects/doc/thesis/figsI-EMIC/', basename];
        system(['mkdir -pv ', dname]);
        system(['cp -v ', fname, ' ', dname, '/.']);
    end
end