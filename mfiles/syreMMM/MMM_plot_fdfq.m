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
figNames{1}  = 'FluxD_3D';
figNames{2}  = 'FluxQ_3D';
figNames{3}  = 'Torque_3D';
figNames{4}  = 'TorRipPP_3D';
figNames{5}  = 'TorRip_3D';
figNames{6}  = 'FluxD_2D';
figNames{7}  = 'FluxQ_2D';
figNames{8}  = 'FluxDQ_2D';
figNames{9}  = 'Torque_2D';
figNames{10} = 'TorRipPP_2D';
figNames{11} = 'TorRip_2D';

for ii=1:length(figNames)
    hfig(ii) = figure();
    figSetting();
    hax(ii) = axes('OuterPosition',[0 0 1 1],...
        'XLim',[min(min(Id)) max(max(Id))],...
        'YLim',[min(min(Iq)) max(max(Iq))],...
        'PlotBoxAspectRatio',[1 1 0.8]);
    xlabel('$i_d$ (A)')
    ylabel('$i_q$ (A)')
    view(3)
    switch ii
        case 1
            zlabel('$\lambda_d$ (Vs)')
            set(gca,'ZLim',[min(min(Fd)) max(max(Fd))])
        case 2
            zlabel('$\lambda_q$ (Vs)')
            set(gca,'ZLim',[min(min(Fq)) max(max(Fq))])
        case 3
            zlabel('$T$ (Nm)')
            set(gca,'ZLim',[min(min(T)) max(max(T))])
        case 4
            zlabel('$\Delta T_{pp}$ (Nm)')
            if (~isnan(max(dTpp(:)))&&~isempty(dTpp)&&max(dTpp(:))~=min(dTpp(:)))
                set(gca,'ZLim',[min(min(dTpp)) max(max(dTpp))])
            else
                set(gca,'ZLim',[0 1]);
            end
        case 5
            zlabel('$\Delta T_{rms}$ (Nm)')
            if (~isnan(max(dT(:)))&&~isempty(dT)&&max(dT(:))~=min(dT(:)))
                set(gca,'ZLim',[min(min(dT)) max(max(dT))])
            else
                set(gca,'ZLim',[0 1]);
            end
        case 6
            view(2)
            title('$\lambda_d$ (Vs)')
        case 7
            view(2)
            title('$\lambda_q$ (Vs)')
        case 8
            view(2)
            title('$|\lambda_{dq}|$ (Vs)')
        case 9
            view(2)
            title('$T$ (Nm)')
        case 10
            view(2)
            title('$\Delta T_{pp}$ (Nm)')
        case 11
            view(2)
            title('$\Delta T_{rms}$ (Nm)')
    end
    set(hfig(ii),'FileName',[pathname resFolder figNames{ii} '.fig'])
    set(hfig(ii),'Name',figNames{ii})
end

surf(hax(1),Id,Iq,Fd,'FaceColor','interp','EdgeColor','none')
contour3(hax(1),Id,Iq,Fd,'EdgeColor','k','ShowText','off')
surf(hax(2),Id,Iq,Fq,'FaceColor','interp','EdgeColor','none')
contour3(hax(2),Id,Iq,Fq,'EdgeColor','k','ShowText','off')
surf(hax(3),Id,Iq,T,'FaceColor','interp','EdgeColor','none')
contour3(hax(3),Id,Iq,T,'EdgeColor','k','ShowText','off')
surf(hax(4),Id,Iq,dTpp,'FaceColor','interp','EdgeColor','none')
contour3(hax(4),Id,Iq,dTpp,'EdgeColor','k','ShowText','off')
surf(hax(5),Id,Iq,dT,'FaceColor','interp','EdgeColor','none')
contour3(hax(5),Id,Iq,dT,'EdgeColor','k','ShowText','off')

contourf(hax(6),Id,Iq,Fd,'ShowText','on');
contourf(hax(7),Id,Iq,Fq,'ShowText','on');
contourf(hax(8),Id,Iq,abs(Fd+j*Fq),'ShowText','on');
contourf(hax(9),Id,Iq,T,'ShowText','on');
contourf(hax(10),Id,Iq,dTpp,'ShowText','on');
contourf(hax(11),Id,Iq,dT,'ShowText','on');

for ii=6:11
    set(hax(ii),'Layer','top','GridColor','k','GridAlpha',1);
end

%% Curves
hfig(12) = figure();
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
xlabel('$i_{dq}$ (A)');
ylabel('$\lambda_{dq}$ (Vs)');
title('Saturation Curves')
set(hfig(12),'FileName',[pathname resFolder 'SaturationCurves.fig'])
set(hfig(12),'Name','SaturationCurves')


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
