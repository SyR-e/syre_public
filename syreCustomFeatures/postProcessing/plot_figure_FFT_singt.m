function [cont,harm] = plot_figure_FFT_singt(x,y,h,xl,yl)

figure()
figSetting(15,12)
subplot(2,1,1)
xlabel(xl)
ylabel(yl)
dy=ceil(max(y)/5);
ymax=max(y);
set(gca,'Xlim',[0 360],'XTick',0:60:360);
ymean=mean(y);
plot(x,y,'-b','LineWidth',2,'DisplayName','waveform')
plot(x,ymean*ones(1,length(x)),'-g','LineWidth',1,'DisplayName',['mean = ' num2str(ymean,3)])
plot(mean(x)*[1 1],[max(y) min(y)],'-rx','LineWidth',1,'DisplayName',['peak-to-peak = ' num2str(max(y)-min(y),3)])
legend('show','Location','South')

subplot(2,1,2)
xlabel('harmonic order')
ylabel([yl ' - harm'])
a=fft(y(1:end-1));
harm=2*abs(a(2:end))/length(y);
cont=abs(a(1))/length(y);
set(gca,'XLim',[-0.5 h+0.5],'XTick',0:1:h);
set(gca,'XLim',[-0.5 h+0.5],'XTick',0:3:48);
bar(1:1:h,[harm(1:h)],'FaceColor','b','BarWidth',0.5)