fontsize = 13;

%% %% PLOT BIFURCATION DIAGRAM ---------------------------------

T2L2no_circ = '/home/erik/Projects/i-emic/test/tuning/solarcont/noO_4deg/EK_Large/0.94+/12h/cdata.txt';
T1L1no_circ = '/home/erik/Projects/i-emic/test/tuning/solarcont/noO_4deg/EK_Large/1.0-/36h/cdata.txt';
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

%% PLOT points P1, P2 in bifurcation diagram

P2circ = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.94+/to0.98/1h/cdata.txt'
P1circ = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.98/flatbot/12h/1h/tdata.txt';

[titles, data] = load_cdata(P1circ);
SIindex = getIndex(titles, 'SIarea');

P1SIarea = data(end,SIindex);
P1SFor   = 0.98;

[titles, data] = load_cdata(P2circ);
SIindex = getIndex(titles, 'SIarea');

P2SIarea = data(end,SIindex);
P2SFor   = 0.98;

hold on;
P1hndl = plot(P1SFor, P1SIarea,'o','color','k','markerfacecolor','k');
P2hndl = plot(P2SFor, P2SIarea,'o','color','k','markerfacecolor','k');
SFline = plot([P1SFor,P2SFor], ylim,'k:');
text(P1SFor, P1SIarea+0.01, ' P_1','horizontalalignment','left','fontsize',fontsize)
text(P2SFor, P2SIarea-0.02, 'P_2 ','horizontalalignment','right','fontsize',fontsize)
hold off;

%% PLOT points L1, L2 in bifurcation diagram


[folds, data] = detectFolds(T2L2circ);


hold on
L2x = folds(1,2);
L2y = data(folds(1,1), SIindex);
L2hndl = plot(L2x, L2y,'o','color','k','markerfacecolor','k');
text(L2x, L2y, '  L_1', 'horizontalalignment','left','fontsize',fontsize);

hold off

%% LEGEND and EXPORTS

legend([no_circ, circ], 'ocean circulation disabled',['ocean ' ...
                    'circulation enabled'], 'location', 'southeast')

xlim([0.9,1.21])
exportfig(['bifdiag4deg12l.eps'], fontsize, [18,26], opts.invert);
!cp -v bifdiag4deg12l.eps /home/erik/Projects/doc/thesis/figsI-EMIC/.

%% Create zoomed plot

% Zoom in on L2
xlim([1.17, 1.185])
ylim([0.907,0.917])
set(legend,'Visible','off');
set(gca,'YaxisLocation','right');

exportfig(['bifdiag4deg12lzL2.eps'], fontsize, [18,12], opts.invert);
!cp -v bifdiag4deg12lzL2.eps /home/erik/Projects/doc/thesis/figsI-EMIC/.

% Zoom in on L1
xlim([0.95, 1.01])
ylim([0.0, 0.20])
set(legend,'Visible','off');
set(gca,'YaxisLocation','right');

exportfig(['bifdiag4deg12lzL1.eps'], fontsize, [18,12], opts.invert);
!cp -v bifdiag4deg12lzL1.eps /home/erik/Projects/doc/thesis/figsI-EMIC/.

axis auto

return

%% %% PLOT FIELDS ---------------------------------

P2SIcirc    = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.94+/to0.98/1h/seaice_output.h5'
P1SIcirc    = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.98/flatbot/12h/1h/seaice_output.h5';

%% Sea ice thickness at P1
opts=[];
opts.invert=false;
opts.exporteps=true;
opts.sifield='H';
plot_seaice(P1SIcirc, opts);
figure(21)
title('');
system('cp -v seaiceH.eps /home/erik/Projects/doc/thesis/figsI-EMIC/P1SIH.eps');

%% Sea ice thickness at P2
opts=[];
opts.invert=false;
opts.exporteps=true;
opts.sifield='H';
plot_seaice(P2SIcirc, opts);
figure(21)
title('');
system('cp -v seaiceH.eps /home/erik/Projects/doc/thesis/figsI-EMIC/P2SIH.eps');


function [SIindex] = getIndex(titles,string)
    SIindex = 0;
    for i = 1:numel(titles)
        if strcmp(titles{i}, string);
            SIindex = i;
            break;
        end
    end
end

function [foldPoints, data] = detectFolds(cdatafile)
    [titles, data] = load_cdata(cdatafile);
    parIndex = getIndex(titles, 'par');
    parData  = data(:,parIndex);
    parDiff  = parData(1:end-1) - parData(2:end);
    sgn = sign(parDiff(1));
    foldPoints = [];
    for i = 1:numel(parDiff)
        if sgn ~= sign(parDiff(i))
            sgn = -sgn;
            foldPoints = [foldPoints ; i, parData(i)];
        end
    end
end
