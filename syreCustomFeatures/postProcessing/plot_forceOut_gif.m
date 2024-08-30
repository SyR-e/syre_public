% Copyright 2024
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_forceOut_gif(forceOut,pathname)

figNames{1} = [pathname 'wf_pressureRadial.gif'];
figNames{2} = [pathname 'wf_pressureTangential.gif'];
figNames{3} = [pathname 'wf_forceRadial.gif'];
figNames{4} = [pathname 'wf_forceTangential.gif'];

for ii=1:length(figNames)
    hfig(ii) = figure();
    figSetting();
    hax(ii) = axes('OuterPosition',[0 0 1 1]);
    switch ii
        case 1
            xlabel('$\xi$ ($^\circ$)')
            ylabel('$p_r$ (Pa)')
            set(gca,'XLim',[0 360],'XTick',0:30:360)
            tmp = max(abs(forceOut.pGapR(:)));
            set(gca,'YLim',tmp*[-1 1]);
            hPlot(ii) = plot(forceOut.thGap,forceOut.pGapR(:,1),'-b');
        case 2
            xlabel('$\xi$ ($^\circ$)')
            ylabel('$p_t$ (Pa)')
            set(gca,'XLim',[0 360],'XTick',0:30:360)
            tmp = max(abs(forceOut.pGapT(:)));
            set(gca,'YLim',tmp*[-1 1]);
            hPlot(ii) = plot(forceOut.thGap,forceOut.pGapT(:,1),'-b');
        case 3
            xlabel('slot')
            ylabel('$F_r$ (N)')
            set(gca,'XLim',[0 length(forceOut.thTooth)+1],'XTick',1:1:length(forceOut.thTooth))
            tmp = max(abs(forceOut.FrTooth(:)),[],'all');
            set(gca,'YLim',tmp*[-1 1]);
            hPlot(ii) = stem(1:1:length(forceOut.thTooth),forceOut.FrTooth(:,1),'-bo','LineWidth',1.5);
        case 4
            xlabel('slot')
            ylabel('$F_t$ (N)')
            set(gca,'XLim',[0 length(forceOut.thTooth)+1],'XTick',1:1:length(forceOut.thTooth))
            tmp = max(abs(forceOut.FtTooth(:)),[],'all');
            set(gca,'YLim',tmp*[-1 1]);
            hPlot(ii) = stem(1:1:length(forceOut.thTooth),forceOut.FtTooth(:,1),'-bo','LineWidth',1.5);
    end
    plot(hax(ii),360*[-1 1],[0 0],'-k','LineWidth',1,'HandleVisibility','off')
    set(hax(ii),'ColorOrderIndex',1)
end

for ii=1:length(forceOut.thRot)
    set(hPlot(1),'YData',forceOut.pGapR(:,ii));
    title(hax(1),['$\theta = ' num2str(forceOut.thRot(ii)) '^\circ$'])
    drawnow
    set(hPlot(2),'YData',forceOut.pGapT(:,ii));
    title(hax(2),['$\theta = ' num2str(forceOut.thRot(ii)) '^\circ$'])
    drawnow
    set(hPlot(3),'YData',forceOut.FrTooth(:,ii));
    title(hax(3),['$\theta = ' num2str(forceOut.thRot(ii)) '^\circ$'])
    drawnow
    set(hPlot(4),'YData',forceOut.FtTooth(:,ii));
    title(hax(4),['$\theta = ' num2str(forceOut.thRot(ii)) '^\circ$'])
    drawnow

    for jj=1:length(hfig)
        % GIF save
        frame = getframe(hfig(jj));
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if ii==1
            imwrite(imind,cm,figNames{jj},'gif','Loopcount',inf,'DelayTime',0.01);
        else
            imwrite(imind,cm,figNames{jj},'gif','WriteMode','append','DelayTime',0.01);
        end
    end
end






