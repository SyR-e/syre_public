function adapt_figure()

% Create figure
figure1 = figure(gcf);
set(figure1,'PaperSize',[20.98 29.68],'Position',[360 400 560 420]);

% Set axes
set(gca,'PlotBoxAspectRatio',[1 0.5 1],...
    'FontWeight','bold');
box('on'); grid('on'); hold('all');

% Set xlabel
xlab = get(gca,'XLabel');
set(xlab,'FontWeight','bold');

% Set ylabel
ylab = get(gca,'YLabel');
set(ylab,'FontWeight','bold');

% Set legend
llab = legend('show');
set(llab,'Location','best');

% Set plot linewidth
plab = get(gca,'Children');

for jj = 1: length(plab)
    set(plab(jj),'LineWidth',2);
end

% set legend
legend1 = legend(gca,'show');
set(legend1,'Location','Best');

% xlim([0 0.5])
% ylim([0 3])

