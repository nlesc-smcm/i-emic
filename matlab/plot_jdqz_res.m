system('grep "residual =" dump | sed "s/.*=//" > residuals');
res = load('residuals');

hold on;
plot(res, 'linewidth',1.5);
ylabel(['$\| \beta \mathbf{A}\mathbf{x} - \alpha \mathbf{B}\mathbf{x} ' ...
        '\| / | \alpha/\beta | $'], 'interpreter','latex')
xlabel('JDQZ step', 'interpreter','latex')
leg = [leg, pwd];
legend(leg);
set(gca,'yscale','log')
hold off;