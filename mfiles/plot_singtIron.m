function plot_singtIron(geo,out,newDir,filemot)

figure
figSetting()
ylabel('[$W$]')

xNames{1} = '$P_{Fe,s,h}$';
xNames{2} = '$P_{Fe,s,e}$';
xNames{3} = '$P_{Fe,r,h}$';
xNames{4} = '$P_{Fe,r,e}$';
xNames{5} = '$P_{PM}$';

set(gca,'XLim',[0.5 5.5],'XTick',1:1:5,'XTickLabel',xNames);

b=bar([out.Pfes_h out.Pfes_c out.Pfer_h out.Pfer_c out.Ppm]);
xtips2 = b(1).XEndPoints;
ytips2 = b(1).YEndPoints;
labels2 = string(round(b(1).YData,2));
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ax = gca;
ylim([0 ax.YTick(end)+ax.YTick(2)])

saveas(gcf,[newDir filemot(1:end-4) '_ironLoss.fig'])
