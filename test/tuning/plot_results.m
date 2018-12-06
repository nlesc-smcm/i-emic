T2L2no_circ = '/home/erik/Projects/i-emic/test/tuning/solarcont/noO_4deg/EK_Large/0.94+/1h/cdata.txt';
T1L1no_circ = '/home/erik/Projects/i-emic/test/tuning/solarcont/noO_4deg/EK_Large/1.0-/1h/cdata.txt';
T2L2circ    = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.94+/36h/cdata.txt';

T1L1circ98p = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.98+/36h/cdata.txt';
T1L1circ98n = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.98-/36h/cdata.txt';

col = lines(5); col(1,:)=[0,0,0];

opts=[];
opts.invert=false;
opts.lsty={'.--','color',col(1,:), 'linewidth',1.2};
opts.plot_entry=15;
opts.point=-1;
opts.hold=false;

[~,~,no_circ] = plot_cdata(T2L2no_circ,    opts)

opts.lsty={'.--','color',col(1,:), 'linewidth',1.2};
opts.hold=true;
plot_cdata(T1L1no_circ, opts)

opts.lsty={'.-','color',col(2,:), 'linewidth',0.5};
opts.hold=true;
[~,~,circ] =plot_cdata(T2L2circ, opts)

opts.lsty={'.-','color',col(2,:), 'linewidth',0.5};
opts.hold=true;
plot_cdata(T1L1circ98p, opts)

opts.lsty={'.-','color',col(2,:), 'linewidth',0.5};
opts.hold=true;
plot_cdata(T1L1circ98n, opts)

title('')
xlabel('Solar Forcing')
ylabel('A^{si} / A')
legend([no_circ, circ], 'ocean circulation disabled',['ocean ' ...
                    'circulation enabled'], 'location', 'southeast')

exportfig(['bifdiag4deg12l.eps'], 12, [18,26], opts.invert);
!cp -v bifdiag4deg12l.eps /home/erik/Projects/doc/thesis/figsI-EMIC/.

% Zoom in on L2
xlim([1.17, 1.185])
ylim([0.907,0.917])
set(legend,'Visible','off');
set(gca,'YaxisLocation','right');

exportfig(['bifdiag4deg12lzL2.eps'], 12, [18,12], opts.invert);
!cp -v bifdiag4deg12lzL2.eps /home/erik/Projects/doc/thesis/figsI-EMIC/.

% Zoom in on L1
xlim([0.95, 1.01])
ylim([0.0, 0.20])
set(legend,'Visible','off');
set(gca,'YaxisLocation','right');

exportfig(['bifdiag4deg12lzL1.eps'], 12, [18,12], opts.invert);
!cp -v bifdiag4deg12lzL1.eps /home/erik/Projects/doc/thesis/figsI-EMIC/.

axis auto