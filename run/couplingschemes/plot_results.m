invert = false;

cfiles = { './0/log_label_decoupled/cdata_lbldecoupled.txt', ...
           './0/log_label_quasicoupled/cdata_lblquasicoupled.txt', ...
           './0/log_label_coupledprecD/cdata_lblcoupledprecD.txt', ...
           './0/log_label_coupledprecB/cdata_lblcoupledprecB.txt', ...
           './0/log_label_coupledprecF/cdata_lblcoupledprecF.txt', ...
           './0/log_label_coupledprecG/cdata_lblcoupledprecG.txt' ...
         };

tfiles = { './0/log_label_time_decoupled/tdata_lbltime_decoupled.txt', ...
           './0/log_label_time_quasicoupled/tdata_lbltime_quasicoupled.txt', ...
           './0/log_label_time_coupledprecB/tdata_lbltime_coupledprecB.txt'
         };
           
lnames = { 'Decoupled', ...
                'Weakly coupled', ...
                'Fully coupled', ...
                'Fully coupled', ...
                'Fully coupled', ...
                'Fully coupled' };

lnamesprec = {  'J: decoupled, M: block diagonal', ...
                'J: weakly coupled, M: block diagonal', ...
                'J: fully coupled, M: block diagonal', ...
                'J: fully coupled, M: backward block GS', ...
                'J: fully coupled, M: forward block GS', ...
                'J: fully coupled, M: forward block GS' };
           
% load continuation data 
N = numel(cfiles);
cdata = cell(N,1);
for i = 1:N
    [c_titles, cdata{i}] = load_cdata(cfiles{i});
end

M = numel(c_titles);
NRfig = 0;
MVfig = 0;
Ffig  = 0;

for i = 2:M
    if strcmp(c_titles{i}, 'NR')
        NRfig = i;
    elseif strcmp(c_titles{i}, 'MV')
        MVfig = i;
    elseif strcmp(c_titles{i}, '||x||')
        Xfig = i;
    elseif strcmp(c_titles{i}, '||F||')
        Ffig = i;
    end    
end

% load transient data
NN = numel(tfiles);
Tdata = cell(NN,1);
for i = 1:NN
    [t_titles, tdata{i}] = load_cdata(tfiles{i});
end

MM = numel(t_titles);
xfig = 0;
dtfig = 0;

for i = 1:MM
    if strcmp(t_titles{i}, 'dt_(y)')
        dtfig = i;
    elseif strcmp(t_titles{i}, '|x|')
        xfig = i;
    elseif strcmp(t_titles{i}, 'NR')
        NRfig = i;
    end    
end

% ---------------------------
% Plot transients ||x||
% ---------------------------
figure(101)
range = 2:NN;
for i = range
    plot(tdata{i}(:,1), tdata{i}(:,xfig),'.-', 'linewidth',0.9)
    hold on;        
end
hold off
grid on
legend(lnames{range})
xlabel('t (y)')
ylabel('||x||_2')
ylim([0,20])
exportfig('transientx.eps', 11,[14,10], invert)
system(['cp -v transientx.eps /home/erik/Projects/doc/thesis/', ...
        'figsNum/.']);

% ---------------------------
% Plot transients dt
% ---------------------------
figure(102)
range = 2:NN;
for i = range
    semilogy(tdata{i}(:,1), tdata{i}(:,dtfig),'.-', 'linewidth',0.9)
    hold on;        
end
hold off
grid on
legend(lnames{range},'location','northeast')
xlabel('t (y)')
ylabel('\Delta t (y)')
ylim([1e-1 5e3])
exportfig('transientdt.eps', 11,[14,10], invert)
system(['cp -v transientdt.eps /home/erik/Projects/doc/thesis/', ...
        'figsNum/.']);

% ---------------------------
% Plot transients NR
% ---------------------------
figure(103)
range = 2:NN;
for i = range
    semilogy(tdata{i}(:,1), tdata{i}(:,NRfig),'.-', 'linewidth',0.9)
    hold on;        
end
hold off
grid on
legend(lnames{range},'location','northeast')
xlabel('t (y)')
ylabel('Newton iterations')
exportfig('transientNR.eps', 11,[14,10], invert)
system(['cp -v transientNR.eps /home/erik/Projects/doc/thesis/', ...
        'figsNum/.']);

% ---------------------------
% Plot newton iterations NR
% ---------------------------
figure(NRfig)
NRrange = 2:3;
for i = NRrange
    plot(cdata{i}(:,1), cdata{i}(:,NRfig),'.-', 'linewidth',0.9)
    hold on;
end
hold off
grid on
legend(lnames{NRrange})
xlabel('Combined forcing')
ylabel('Newton iterations')
exportfig('NRcomparison.eps', 11,[14,10], invert)
system(['cp -v NRcomparison.eps /home/erik/Projects/doc/thesis/', ...
        'figsNum/.']);

% ---------------------------
% Plot x norm
% ---------------------------
figure(Xfig)
xrange = 2:3;
for i = xrange
    p = plot(cdata{i}(:,1), cdata{i}(:,Xfig),'.-', 'linewidth',0.9)
    hold on;
end
hold off
grid on
legend(lnames{xrange},'location','northwest')
uistack(p,'down',1)
xlabel('Combined forcing')
ylabel('||x||_2')
exportfig('xcomparison.eps', 11,[14,10], invert)
system(['cp -v xcomparison.eps /home/erik/Projects/doc/thesis/', ...
        'figsNum/.']);

% ---------------------------
% Plot residual
% ---------------------------
figure(Ffig)
Frange = 2:3;
for i = Frange
    semilogy(cdata{i}(:,1), cdata{i}(:,Ffig),'.-', 'linewidth',0.9)
    hold on;
end
hold off
grid on
legend(lnames{Frange},'location','northeast')
xlabel('Combined forcing')
ylabel('||F||_2')
exportfig('Fcomparison.eps', 11,[14,10], invert)
system(['cp -v Fcomparison.eps /home/erik/Projects/doc/thesis/', ...
        'figsNum/.']);

% --------------------------------
% Plot matrix vector products
% --------------------------------
figure(MVfig)
MVrange = 2:N-1;
for i = MVrange
    plot(cdata{i}(:,1), cdata{i}(:,MVfig),'.-', 'linewidth',0.9)
    hold on;
end
hold off
grid on
legend(lnamesprec{MVrange},'location','northwest')
xlabel('Combined forcing')
ylabel('FGMRES iterations')
ylim([0,30]);
exportfig('MVcomparison.eps', 11,[14,10], invert)
system(['cp -v MVcomparison.eps /home/erik/Projects/doc/thesis/', ...
        'figsNum/.']);

return












[c_titles, cdata] = plot_cdata([], 'r--');

[c_titles, cdata] = plot_cdata(['./0/log_label_testB/' ...
                    'cdata_lbltestB.txt'], 'b-');

[c_titles, cdata] = plot_cdata(['./0/log_label_testC2/' ...
                    'cdata_lbltestC2.txt'], 'r-');

[c_titles, cdata] = plot_cdata(['./0/log_label_testC3/' ...
                    'cdata_lbltestC3.txt'], 'k-');


for i = 2:N
    figure(i);
    hold off;
    if (i == MVfig)
        set(gca, 'yscale', 'log')
    end
    legend('Decoupled', 'Quasi-coupled', 'Fully coupled  (prec: B-GS)', ...
           'Fully coupled  (prec: D)','Fully coupled  (prec: F-GS)', ...
           'location','southeast')
    ylabel(c_titles{i})
    title('')

end

