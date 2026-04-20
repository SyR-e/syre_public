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
if strcmp(motorModel.data.motorType,'EE')
    Ir = motorModel.FluxMap_dq.Ir;
    Fr = motorModel.FluxMap_dq.Fr;
end

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = checkPathSyntax([motName '_results\MMM results\' 'dq Flux Maps - ' int2str(motorModel.data.tempPM) 'deg\']);

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
        'XLim',[min(min(min(Id))) max(max(max(Id)))],...
        'YLim',[min(min(min(Iq))) max(max(max(Iq)))],...
        'PlotBoxAspectRatio',[1 1 0.8]);
    xlabel('$i_d$ (A)')
    ylabel('$i_q$ (A)')
    view(3)
    switch ii
        case 1
            zlabel('$\lambda_d$ (Vs)')
            set(gca,'ZLim',[min(Fd(:)) max(Fd(:))])
        case 2
            zlabel('$\lambda_q$ (Vs)')
            set(gca,'ZLim',[min(Fq(:)) max(Fq(:))])
        case 3
            zlabel('$T$ (Nm)')
            set(gca,'ZLim',[min(T(:)) max(T(:))])
        case 4
            zlabel('$\Delta T_{pp}$ (Nm)')
            if (~isnan(max(dTpp(:)))&&~isempty(dTpp)&&max(dTpp(:))~=min(dTpp(:)))
                set(gca,'ZLim',[min(dTpp(:)) max(dTpp(:))])
            else
                set(gca,'ZLim',[0 1]);
            end
        case 5
            zlabel('$\Delta T_{rms}$ (Nm)')
            if (~isnan(max(dT(:)))&&~isempty(dT)&&max(dT(:))~=min(dT(:)))
                set(gca,'ZLim',[min(dT(:)) max(dT(:))])
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

if ~strcmp(motorModel.data.motorType,'EE')
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
    h_count = 11;
else
    surf(hax(1),Id(:,:,end),Iq(:,:,end),Fd(:,:,end),'FaceColor','interp','EdgeColor','none')
    contour3(hax(1),Id(:,:,end),Iq(:,:,end),Fd(:,:,end),'EdgeColor','k','ShowText','off')
    surf(hax(2),Id(:,:,end),Iq(:,:,end),Fq(:,:,end),'FaceColor','interp','EdgeColor','none')
    contour3(hax(2),Id(:,:,end),Iq(:,:,end),Fq(:,:,end),'EdgeColor','k','ShowText','off')
    surf(hax(3),Id(:,:,end),Iq(:,:,end),T(:,:,end),'FaceColor','interp','EdgeColor','none')
    contour3(hax(3),Id(:,:,end),Iq(:,:,end),T(:,:,end),'EdgeColor','k','ShowText','off')
    surf(hax(4),Id(:,:,end),Iq(:,:,end),dTpp(:,:,end),'FaceColor','interp','EdgeColor','none')
    contour3(hax(4),Id(:,:,end),Iq(:,:,end),dTpp(:,:,end),'EdgeColor','k','ShowText','off')
    surf(hax(5),Id(:,:,end),Iq(:,:,end),dT(:,:,end),'FaceColor','interp','EdgeColor','none')
    contour3(hax(5),Id(:,:,end),Iq(:,:,end),dT(:,:,end),'EdgeColor','k','ShowText','off')
    
    contourf(hax(6),Id(:,:,end),Iq(:,:,end),Fd(:,:,end),'ShowText','on');
    contourf(hax(7),Id(:,:,end),Iq(:,:,end),Fq(:,:,end),'ShowText','on');
    contourf(hax(8),Id(:,:,end),Iq(:,:,end),abs(Fd(:,:,end)+j*Fq(:,:,end)),'ShowText','on');
    contourf(hax(9),Id(:,:,end),Iq(:,:,end),T(:,:,end),'ShowText','on');
    contourf(hax(10),Id(:,:,end),Iq(:,:,end),dTpp(:,:,end),'ShowText','on');
    contourf(hax(11),Id(:,:,end),Iq(:,:,end),dT(:,:,end),'ShowText','on');
    % Multi Ir plot
    size_Ir = size(Ir(1,1,:));
    Ir_plot_index = round(linspace(1, size_Ir(end), 10)); % 10 è il # di plot tra da vedere
    Ir_plot_index(1) = 1;
    Ir_plot_index(end) = max(size_Ir);
    Ir_plot_index = unique(Ir_plot_index); %Eventuali indici ripetuti vengono rimossi
    Ir_values = Ir(1,1,:);
    Ir_plot_values = Ir_values(Ir_plot_index);
    hfig(12) = surfplot_slider(Id(:,:,end),Iq(:,:,end),T,Ir_plot_values,[pathname resFolder 'Torque_Ir.fig']);
    h_count = 12;
end

for ii=6:11
    set(hax(ii),'Layer','top','GridColor','k','GridAlpha',1);
end

%% Curves
if ~strcmp(motorModel.data.motorType,'EE')
    hfig(h_count+1) = figure();
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
    set(hfig(h_count+1),'FileName',[pathname resFolder 'SaturationCurves.fig'])
    set(hfig(h_count+1),'Name','SaturationCurves')
else
    hfig(h_count+1) = figure();
    figSetting()
    [~, index] = min(abs(Iq(:,1,end)));
    plot(Id(index,:,end),Fd(index,:,end),'Color',[0 0 0.9],'LineStyle','-','DisplayName','$\lambda_d(i_d,0)$')
    [~, index] = min(abs(Id(1,:,end)));
    plot(Iq(:,index,end),Fq(:,index,end),'Color',[0 0.9 0],'LineStyle','-','DisplayName','$\lambda_q(0,i_q)$')
    plot(Id(1,:,end),Fd(1,:,end),'Color',[0 0 0.9],'LineStyle',':','DisplayName','$\lambda_d(i_d,I_{q,min})$'),
    plot(Id(end,:,end),Fd(end,:,end),'Color',[0 0 0.9],'LineStyle','--','DisplayName','$\lambda_d(i_d,I_{q,max})$')
    plot(Iq(:,1,end),Fq(:,1,end),'Color',[0 0.9 0],'LineStyle',':','DisplayName','$\lambda_q(I_{d,min},i_q)$')
    plot(Iq(:,end,end),Fq(:,end,end),'Color',[0 0.9 0],'LineStyle','--','DisplayName','$\lambda_q(I_{d,max},i_q)$')
    legend('show','Location','NorthWest');
    xlabel('$i_{dq}$ (A)');
    ylabel('$\lambda_{dq}$ (Vs)');
    title('Saturation Curves')
    set(hfig(h_count+1),'FileName',[pathname resFolder 'SaturationCurves.fig'])
    set(hfig(h_count+1),'Name','SaturationCurves')
end


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
