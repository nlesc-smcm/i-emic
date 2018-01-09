files = { './0/log_label_decoupled/cdata_lbldecoupled.txt', ...
          './0/log_label_quasicoupled/cdata_lblquasicoupled.txt', ...
          './0/log_label_coupledprecD/cdata_lblcoupledprecD.txt', ...
          './0/log_label_coupledprecB/cdata_lblcoupledprecB.txt', ...
          './0/log_label_coupledprecF/cdata_lblcoupledprecF.txt', ...
          './0/log_label_coupledprecG/cdata_lblcoupledprecG.txt' };

lnames = { 'Decoupled', ...
                'Quasi-coupled', ...
                'Coupled', ...
                'Coupled', ...
                'Coupled', ...
                'Coupled' };

lnamesprec = { 'Decoupled', ...
                'Quasi-coupled', ...
                'Coupled: prec D', ...
                'Coupled: prec B', ...
                'Coupled: prec F', ...
                'Coupled: prec G' };

invert = true;
N = numel(files);
cdata  = cell(N,1);
for i = 1:N
    [titles, cdata{i}] = load_cdata(files{i});
end

M = numel(titles);
NRfig = 0;
MVfig = 0;
Ffig  = 0;

for i = 2:M
    if strcmp(titles{i}, 'NR')
        NRfig = i;
    elseif strcmp(titles{i}, 'MV')
        MVfig = i;
    elseif strcmp(titles{i}, '||F||')
        Ffig = i;
    end    
end

% ---------------------------
% Plot newton iterations NR
% ---------------------------
figure(NRfig)
NRrange = 1:3;
for i = NRrange
    plot(cdata{i}(:,1), cdata{i}(:,NRfig),'.-', 'linewidth',1.3)
    hold on;
end
hold off
grid on
legend(lnames{NRrange})
xlabel('Continuation parameter')
ylabel('Newton iterations')
exportfig('NRcomparison.eps', 12, [14,10], invert)

% ---------------------------
% Plot residual
% ---------------------------
figure(Ffig)
Frange = 1:3;
for i = Frange
    semilogy(cdata{i}(:,1), cdata{i}(:,Ffig),'.-', 'linewidth',1.3)
    hold on;
end
hold off
grid on
legend(lnames{Frange})
xlabel('Continuation parameter')
ylabel('||F||_2')
exportfig('Fcomparison.eps', 12, [14,10], invert)

% --------------------------------
% Plot matrix vector products
% --------------------------------
figure(MVfig)
MVrange = 1:N;
for i = MVrange
    plot(cdata{i}(:,1), cdata{i}(:,MVfig),'.-', 'linewidth',1.3)
    hold on;
end
hold off
grid on
legend(lnamesprec{MVrange},'location','southeast')
xlabel('Continuation parameter')
ylabel('FGMRES iterations')
exportfig('MVcomparison.eps', 12, [14,10], invert)

return












[titles, cdata] = plot_cdata([], 'r--');

[titles, cdata] = plot_cdata(['./0/log_label_testB/' ...
                    'cdata_lbltestB.txt'], 'b-');

[titles, cdata] = plot_cdata(['./0/log_label_testC2/' ...
                    'cdata_lbltestC2.txt'], 'r-');

[titles, cdata] = plot_cdata(['./0/log_label_testC3/' ...
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
    ylabel(titles{i})
    title('')

end

