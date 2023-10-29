% Copyright 2020
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

function MMM_plot_fdfq(motorModel)

% load data
Id   = motorModel.FluxMap_dq.Id;
Iq   = motorModel.FluxMap_dq.Iq;
Fd   = motorModel.FluxMap_dq.Fd;
Fq   = motorModel.FluxMap_dq.Fq;
T    = motorModel.FluxMap_dq.T;
dT   = motorModel.FluxMap_dq.dT;
dTpp = motorModel.FluxMap_dq.dTpp;

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'dq Flux Maps - ' int2str(motorModel.data.tempPM) 'deg\'];

%% Surfaces
figNames{1} = 'FluxD';
figNames{2} = 'FluxQ';
figNames{3} = 'Torque';
figNames{4} = 'TorRipPP';
figNames{5} = 'TorRip';

for ii=1:length(figNames)
    hfig(ii) = figure();
    figSetting();
    hax(ii) = axes('OuterPosition',[0 0 1 1],...
        'XLim',[min(min(Id)) max(max(Id))],...
        'YLim',[min(min(Iq)) max(max(Iq))],...
        'PlotBoxAspectRatio',[1 1 0.8]);
    xlabel('$i_d$ [A]')
    ylabel('$i_q$ [A]')
    view(3)
    switch ii
        case 1
            zlabel('$\lambda_d$ [Vs]')
            set(gca,'ZLim',[min(min(Fd)) max(max(Fd))])
        case 2
            zlabel('$\lambda_q$ [Vs]')
            set(gca,'ZLim',[min(min(Fq)) max(max(Fq))])
        case 3
            zlabel('$T$ [Nm]')
            set(gca,'ZLim',[min(min(T)) max(max(T))])
        case 4
            zlabel('$\Delta T_{pp}$ [Nm]')
            if (~isnan(max(dTpp(:)))&&~isempty(dTpp)&&max(dTpp(:))~=min(dTpp(:)))
                set(gca,'ZLim',[min(min(dTpp)) max(max(dTpp))])
            else
                set(gca,'ZLim',[0 1]);
            end
        case 5
            zlabel('$\Delta T_{rms}$ [Nm]')
            if (~isnan(max(dT(:)))&&~isempty(dT)&&max(dT(:))~=min(dT(:)))
                set(gca,'ZLim',[min(min(dT)) max(max(dT))])
            else
                set(gca,'ZLim',[0 1]);
            end
    end
    set(hfig(ii),'FileName',[pathname resFolder figNames{ii} '.fig'])
    set(hfig(ii),'Name',figNames{ii})
end

surf(hax(1),Id,Iq,Fd,'FaceColor','interp','EdgeColor','interp')
surf(hax(2),Id,Iq,Fq,'FaceColor','interp','EdgeColor','interp')
surf(hax(3),Id,Iq,T,'FaceColor','interp','EdgeColor','interp')
surf(hax(4),Id,Iq,dTpp,'FaceColor','interp','EdgeColor','interp')
surf(hax(5),Id,Iq,dT,'FaceColor','interp','EdgeColor','interp')

%% Curves
hfig(6) = figure();
figSetting()
[~, index] = min(abs(Iq(:,1)));
plot(Id(index,:),Fd(index,:),'Color',[0 0 0.9],'LineStyle','-','DisplayName','$\lambda_d(i_d,0)$')
[~, index] = min(abs(Id(1,:)));
plot(Iq(:,index),Fq(:,index),'Color',[0 0.9 0],'LineStyle','-','DisplayName','$\lambda_q(0,i_q)$')
plot(Id(1,:),Fd(1,:),'Color',[0 0 0.9],'LineStyle',':','DisplayName','$\lambda_d(i_d,I_{q,min})$'),
plot(Id(end,:),Fd(end,:),'Color',[0 0 0.9],'LineStyle','--','DisplayName','$\lambda_d(i_d,I_{q,max})$')
plot(Iq(:,1),Fq(:,1),'Color',[0 0.9 0],'LineStyle',':','DisplayName','$\lambda_q(I_{d,min},i_q)$')
plot(Iq(:,end),Fq(:,end),'Color',[0 0.9 0],'LineStyle','--','DisplayName','$\lambda_q(I_{d,max},i_q)$')
legend('show','Location','NorthWest');
xlabel('$i_{dq}$ [A]');
ylabel('$\lambda_{dq}$ [Vs]');
title('Magnetic Model')
set(hfig(6),'FileName',[pathname resFolder 'Curves.fig'])



%% Save figures
answer = 'No';
answer = questdlg('Save figures?','Save','Yes','No',answer);
if strcmp(answer,'Yes')
    if ~exist([pathname resFolder],'dir')
        mkdir([pathname resFolder]);
    end
    
    for ii=1:length(hfig)
        savePrintFigure(hfig(ii));
    end
end
