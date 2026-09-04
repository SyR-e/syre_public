function ELUkg = plot_ELUkg_weights()


[~,~,ELUkg] = calc_ELU();

names = fieldnames(ELUkg);

values = [];

for ii=1:length(names)
    values(ii) = ELUkg.(names{ii});
end

xNames{1} = 'FeSi';
xNames{2} = 'FeCo';
xNames{3} = 'Cu';
xNames{4} = 'Al';
xNames{5} = 'NdFeB';
xNames{6} = 'Dy-free NdFeB';
xNames{7} = 'SmCo';
xNames{8} = 'Ferrite';

colors.lamination = [0.5 0.5 0.5];
colors.conductors = [1.0 0.5 0.0];
colors.magnets    = [0.0 0.0 1.0];
colors.sleeve     = [0.0 0.6 0.0];

figure();
figSetting(14,10)
ylabel('ELU/kg')

set(gca,'XLim',[0.5 8.5],'XTick',1:1:9,'XTickLabel',xNames)
set(gca,'YLim',[0 400],'YTick',0:50:400);

X = [0.5 2.5 2.5 0.5];
Y = [0 0 1000 1000];

hFill = fill(X,Y,'r');
set(hFill,'FaceColor',colors.lamination,'FaceAlpha',0.3,'EdgeColor','none')

X = [2.5 4.5 4.5 2.5];

hFill = fill(X,Y,'r');
set(hFill,'FaceColor',colors.conductors,'FaceAlpha',0.3,'EdgeColor','none')

X = [4.5 8.5 8.5 4.5];

hFill = fill(X,Y,'r');
set(hFill,'FaceColor',colors.magnets,'FaceAlpha',0.3,'EdgeColor','none')


hBar = bar(values);
set(hBar,'FaceColor',[1 0 0],'EdgeColor',0.5*[1 0 0],'BarWidth',0.5)

text(1.5,360,'Lamination','VerticalAlignment','bottom','HorizontalAlignment','center','Color',colors.lamination*0.5)
text(3.5,360,'Conductor','VerticalAlignment','bottom','HorizontalAlignment','center','Color',colors.conductors*0.5)
text(6.5,360,'Magnet','VerticalAlignment','bottom','HorizontalAlignment','center','Color',colors.magnets*0.5)


for ii=1:length(names)
    text(ii,values(ii),num2str(round(values(ii),1)),'VerticalAlignment','bottom','HorizontalAlignment','center','Color',0.5*[1 0 0])
end



