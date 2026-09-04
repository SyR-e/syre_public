function [cont,harm] = plot_figure_FFT_singt(x,y,h,xl,yl)

figure()
figSetting(15,12)
subplot(2,1,1)
xlabel(xl)
ylabel(yl)
dy=ceil(max(y)/5);
ymax=max(y);
set(gca,'Xlim',[0 360],'XTick',0:60:360);
yLim = [(min(y)) (max(y))]*1.1;
if sign(max(yLim))==sign(min(yLim))
    if sign(yLim)==+1
        yLim(1) = 0;
    else
        yLim(2) = 0;
    end
end
ylim(yLim)
ymean=mean(y);
plot(x,y,'-','DisplayName','waveform')
plot(x,ymean*ones(1,length(x)),':','DisplayName',['mean = ' num2str(ymean,3)])
plot(mean(x)*[1 1],[max(y) min(y)],':x','DisplayName',['peak-to-peak = ' num2str(max(y)-min(y),3)])
legend('show','Location','southeast')

subplot(2,1,2)
xlabel('harmonic order')
ylabel([yl ' - harm'])
a=fft(y(1:end-1));
harm=2*abs(a(2:end))/length(y);
cont=abs(a(1))/length(y);
set(gca,'XLim',[-0.5 h+0.5],'XTick',0:6:48);
set(gca,'YLim',[0 max(harm)*1.1]);
bar(1:1:h,[harm(1:h)],'FaceColor','b','BarWidth',1)