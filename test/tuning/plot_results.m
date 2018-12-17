fontsize = 13;

%% %% PLOT BIFURCATION DIAGRAM ---------------------------------

T2L2no_circ = '/home/erik/Projects/i-emic/test/tuning/solarcont/noO_4deg/EK_Large/0.94+/12h/cdata.txt';
T1L1no_circ = '/home/erik/Projects/i-emic/test/tuning/solarcont/noO_4deg/EK_Large/1.0-/36h/cdata.txt';
T2L2circ    = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.94+/36h/cdata.txt';

T1L1circ98p = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.98+/36h/cdata.txt';
T1L1circ98n = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.98-/5d/cdata.txt';

T1L1circ97468n = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.97468-/1h4/cdata.txt';

T1L1circ975n = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.975-/5d/cdata.txt';
T1L1circ975p = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.975+/5d/cdata.txt';

T1L1no_circ98n  = '/home/erik/Projects/i-emic/test/tuning/solarcont/noO_4deg/EK_Large/0.98-/5d/cdata.txt';
T1L1no_circ98n2 = '/home/erik/Projects/i-emic/test/tuning/solarcont/noO_4deg/EK_Large/0.98-/5d/5d/cdata.txt';
T1L1no_circ98p  = '/home/erik/Projects/i-emic/test/tuning/solarcont/noO_4deg/EK_Large/0.98+/5d/cdata.txt';
T1L1no_circ98p2 = '/home/erik/Projects/i-emic/test/tuning/solarcont/noO_4deg/EK_Large/0.98+/5d/5d/cdata.txt';

L1no_circ97637n = '/home/erik/Projects/i-emic/test/tuning/solarcont/noO_4deg/EK_Large/0.97637-/ev/1h/cdata.txt';

col = lines(5); col(1,:)=[0,0,0];

lsty_no_circ = '.-';
lsty_circ    = '.-';
lw_circ = 0.8;
lw_no_circ = 0.8;

opts=[];
opts.invert=false;
opts.lsty={lsty_no_circ,'color',col(1,:), 'linewidth',lw_no_circ};
opts.plot_entry=15;
opts.point=-1;
opts.hold=false;

[~,~,no_circ] = plot_cdata(T2L2no_circ, opts);

opts.lsty={lsty_no_circ,'color',col(1,:), 'linewidth',lw_no_circ};
opts.hold=true;
plot_cdata(T1L1no_circ, opts);

opts.lsty={lsty_no_circ,'color',col(1,:), 'linewidth',lw_no_circ};
opts.hold=true;
plot_cdata(T1L1no_circ98n, opts);

opts.lsty={lsty_no_circ,'color',col(1,:), 'linewidth',lw_no_circ};
opts.hold=true;
plot_cdata(T1L1no_circ98n2, opts);

opts.lsty={lsty_no_circ,'color',col(1,:), 'linewidth',lw_no_circ};
opts.hold=true;
plot_cdata(T1L1no_circ98p, opts);

opts.lsty={lsty_no_circ,'color',col(1,:), 'linewidth',lw_no_circ};
opts.hold=true;
plot_cdata(T1L1no_circ98p2, opts);

opts.lsty={lsty_circ,'color',col(2,:), 'linewidth',lw_circ};
opts.hold=true;
[~,~,circ] =plot_cdata(T2L2circ, opts);

opts.lsty={lsty_circ,'color',col(2,:), 'linewidth',lw_circ};
opts.hold=true;
plot_cdata(T1L1circ98p, opts);

opts.lsty={lsty_circ,'color',col(2,:), 'linewidth',lw_circ};
opts.hold=true;
plot_cdata(T1L1circ98n, opts);

opts.lsty={lsty_circ,'color',col(2,:), 'linewidth',lw_circ};
opts.hold=true;
plot_cdata(T1L1circ97468n, opts);

opts.lsty={lsty_circ,'color',col(2,:), 'linewidth',lw_circ};
opts.hold=true;
plot_cdata(T1L1circ975n, opts);

opts.lsty={lsty_circ,'color',col(2,:), 'linewidth',lw_circ};
opts.hold=true;
plot_cdata(T1L1circ975p, opts);

title('')
xlabel('Solar forcing \lambda_\Sigma')
ylabel('A^{si} / A')

%% PLOT points P1, P2 in bifurcation diagram

P2circ    = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.94+/to0.98/1h/cdata.txt'
P1circ    = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.98/flatbot/12h/1h/tdata.txt';
P1no_circ = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.98/no_ocean/6h/tdata.txt';


[titles, data]  =  load_cdata(P1circ);
SIindex         =  getIndex(titles, 'SIarea');
P1SIarea_circ   =  data(end,SIindex);
P1SFor          =  0.98;

[titles, data]   = load_cdata(P1no_circ);
SIindex          = getIndex(titles, 'SIarea');
P1SIarea_no_circ = data(end,SIindex);
P1SFor           = 0.98;

[titles, data] = load_cdata(P2circ);
SIindex        = getIndex(titles, 'SIarea');
P2SIarea       = data(end,SIindex);
P2SFor         = 0.98;

hold on;
P1hndla = plot(P1SFor, P1SIarea_circ,'o','color','k','markerfacecolor','k');
P1hndlb = plot(P1SFor, P1SIarea_no_circ,'o','color','k','markerfacecolor','k');
P2hndl = plot(P2SFor, P2SIarea,'o','color','k','markerfacecolor','k');
SFline = plot([P1SFor,P2SFor], ylim,'k:');
text(P1SFor, P1SIarea_no_circ+0.005, ' P_1^a','horizontalalignment','left','fontsize',fontsize)
text(P1SFor, P1SIarea_circ+0.005, ' P_1^b','horizontalalignment','left','fontsize',fontsize)
text(P2SFor, P2SIarea-0.02, 'P_2 ','horizontalalignment','right','fontsize',fontsize)
hold off;

%% PLOT points L1, L2 in bifurcation diagram

hold on

[folds, data] = detectFolds(L1no_circ97637n);
L1xa = folds(1,2)
L1ya = data(folds(1,1), SIindex)
L1hndla = plot(L1xa, L1ya,'o','color','k','markerfacecolor','w');
text(L1xa, L1ya, 'L^a_1 ', 'horizontalalignment','right','fontsize',fontsize);

[folds, data] = detectFolds(T1L1circ975n);
idx = find(folds(:,2) == min(folds(:,2)));
L1xb = folds(idx,2);
L1yb = data(folds(idx,1), SIindex);
L1hndlb = plot(L1xb, L1yb,'o','color','k','markerfacecolor','w');
text(L1xb, L1yb, 'L^b_1 ', 'horizontalalignment','right','fontsize',fontsize);

[folds, data] = detectFolds(T2L2circ);
L1xb = min(folds(:,2));


L2x = folds(1,2);
L2y = data(folds(1,1), SIindex);
L2hndl = plot(L2x, L2y,'o','color','k','markerfacecolor','w');
text(L2x, L2y, '  L_2', 'horizontalalignment','left','fontsize',fontsize);

hold off
%% LEGEND and EXPORTS

legend([no_circ, circ], 'ocean circulation disabled (a)',['ocean ' ...
                    'circulation enabled (b)'], 'location', 'southeast')


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
xlim([0.972, 0.988])
ylim([0.01, 0.08])
set(legend,'Visible','off');
set(gca,'YaxisLocation','right');

exportfig(['bifdiag4deg12lzL1.eps'], fontsize, [18,12], opts.invert);
!cp -v bifdiag4deg12lzL1.eps /home/erik/Projects/doc/thesis/figsI-EMIC/.

axis auto

%% %% PLOT FIELDS ---------------------------------

maskFile    = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.98/flatbot/12h/fort.44';
rstrFile    = '/home/erik/Projects/i-emic/data/mkmask/mask.glo_atl';

P1SIcirc    = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.98/flatbot/12h/1h/seaice_output.h5';
P1SIno_circ = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.98/no_ocean/6h/seaice_output.h5';

P1OCcirc    = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.98/flatbot/12h/ocean_output.h5';
P1OCno_circ = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.98/no_ocean/6h/ocean_output.h5';

P1ATcirc    = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.98/flatbot/12h/atmos_output.h5';
P1ATno_circ = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.98/no_ocean/6h/atmos_output.h5';

P2SIcirc    = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.94+/to0.98/1h/seaice_output.h5';
P2OCcirc    = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.94+/to0.98/1h/ocean_output.h5';
P2ATcirc    = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.94+/to0.98/1h/atmos_output.h5';

L1OCcirc    = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.9747/flatbot/12h/ocean_output.h5';
L1ATcirc    = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.9747/flatbot/12h/atmos_output.h5';
L1SIcirc    = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.9747/flatbot/12h/seaice_output.h5';

L1OCno_circ = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.976373/no_ocean/12h/ocean_output.h5';
L1ATno_circ = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.976373/no_ocean/12h/atmos_output.h5';
L1SIno_circ = '/home/erik/Projects/i-emic/test/tuning/transient/4deg12layers/sol0.976373/no_ocean/12h/seaice_output.h5';

L2OCcirc = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.94+/evs_L2/shift100/ocean_output.h5';
L2ATcirc = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.94+/evs_L2/shift100/atmos_output.h5';
L2SIcirc = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.94+/evs_L2/shift100/seaice_output.h5';

%% Fields at P1a and P1b

opts.mstream=false;
opts.sst=true;
opts.Ta=true;
opts.maskfile = maskFile;
opts.invert=false;
opts.exporteps=false;
opts.sifield='H';
opts.fig_ctr=1;
opts.only_contour=true;

% ----------------------
SIFiles  = {'P1bSI.eps', 'P1aSI.eps', 'P2SI.eps', 'L1aSI.eps', ...
            'L1bSI.eps', 'L2SI.eps'};
pointsSI = { P1SIcirc, P1SIno_circ, P2SIcirc, L1SIno_circ, L1SIcirc, ...
             L2SIcirc};
MOCFiles = {'P1bMOC.eps', 'P1aMOC.eps', 'P2MOC.eps', 'L1aMOC.eps', ...
            'L1bMOC.eps', 'L2MOC.eps'};

AMOCFiles = {'P1bAMOC.eps', 'P1aAMOC.eps', 'P2AMOC.eps', 'L1aAMOC.eps', ...
            'L1bAMOC.eps', 'L2AMOC.eps'};

pointsOC = { P1OCcirc, P1OCno_circ, P2OCcirc, L1OCno_circ, L1OCcirc, ...
             L2OCcirc};

TaFiles = {'P1bTa.eps', 'P1aTa.eps', 'P2Ta.eps', 'L1aTa.eps', ...
            'L1bTa.eps', 'L2Ta.eps'};

pointsAT = { P1ATcirc, P1ATno_circ, P2ATcirc, L1ATno_circ, L1ATcirc, ...
             L2ATcirc};

for i = 1:numel(pointsSI)
    plot_seaice(pointsSI{i}, opts);
    figure(1); hold on;
    plot_ocean(pointsOC{i}, opts);
    hold off;

    title('');
    exportfig(SIFiles{i}, fontsize, [20,10])
    system(['cp -v ', SIFiles{i}, ' /home/erik/Projects/doc/thesis/figsI-EMIC/.']);
end

for i = 1:numel(pointsOC)
    opts.restrict_sol=false;        
    opts.sst=false;
    opts.mstream=true;
    plot_ocean(pointsOC{i}, opts);
    title('');
    exportfig(MOCFiles{i}, fontsize, [20,10])
    system(['cp -v ', MOCFiles{i}, ' /home/erik/Projects/doc/thesis/figsI-EMIC/.']);

    opts.restrict_sol=true;
    opts.rmask_file=rstrFile;
    plot_ocean(pointsOC{i}, opts);
    title('');
    exportfig(AMOCFiles{i}, fontsize, [20,10])
    system(['cp -v ', AMOCFiles{i}, ' /home/erik/Projects/doc/thesis/figsI-EMIC/.']);
end

for i = 1:numel(pointsAT)
    plot_atmos(pointsAT{i}, opts);
    title('')
    exportfig(TaFiles{i}, fontsize, [20,10])
    system(['cp -v ', TaFiles{i}, ' /home/erik/Projects/doc/thesis/figsI-EMIC/.']);    
end


%% Eigenvectors at L1a
L1ATno_circ = '/home/erik/Projects/i-emic/test/tuning/solarcont/noO_4deg/EK_Large/0.97637-/ev/1h/ev_step_4.1.h5';
L1SIno_circ = '/home/erik/Projects/i-emic/test/tuning/solarcont/noO_4deg/EK_Large/0.97637-/ev/1h/ev_step_4.2.h5';
opts.readEV = true;
opts.Ta = true;
opts.evindex=0;
plot_atmos(L1ATno_circ,opts);
title('')
exportfig('L1aTa_ev0.eps',fontsize,[20,10])
system('cp -v L1aTa_ev0.eps /home/erik/Projects/doc/thesis/figsI-EMIC/L1aTa_ev0.eps');

opts.evindex=1;
plot_atmos(L1ATno_circ,opts);
title('')
exportfig('L1aTa_ev1.eps',fontsize,[20,10])
system('cp -v L1aTa_ev1.eps /home/erik/Projects/doc/thesis/figsI-EMIC/L1aTa_ev1.eps');

%% Eigenvectors at L2
L2OCcirc = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.94+/evs_L2/shift100/ev_step_10.0.h5';
L2ATcirc = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.94+/evs_L2/shift100/ev_step_10.1.h5';
L2SIcirc = '/home/erik/Projects/i-emic/test/tuning/solarcont/full_4deg/EK_Large/0.94+/evs_L2/shift100/ev_step_10.2.h5';

opts.evindex=0;
plot_atmos(L2ATcirc,opts);
title('')
exportfig('L2Ta_ev0.eps', fontsize, [20,10])
system('cp -v L2Ta_ev0.eps /home/erik/Projects/doc/thesis/figsI-EMIC/L2Ta_ev0.eps');


%% FUNCTIONS ----------------------------------------
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
