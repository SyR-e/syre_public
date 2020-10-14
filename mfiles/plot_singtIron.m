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

bar([out.Pfes_h out.Pfes_c out.Pfer_h out.Pfer_c out.Ppm]);

saveas(gcf,[newDir filemot(1:end-4) '_ironLoss.fig'])
