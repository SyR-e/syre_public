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

function MMM_plot_idiq(motorModel)

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = checkPathSyntax([motName '_results\MMM results\' 'Inverse Model dq - ' int2str(motorModel.data.tempPM) 'deg\']);

% Load data
Id = motorModel.FluxMapInv_dq.Id;
Iq = motorModel.FluxMapInv_dq.Iq;
Fd = motorModel.FluxMapInv_dq.Fd;
Fq = motorModel.FluxMapInv_dq.Fq;
T  = motorModel.FluxMapInv_dq.T;
if strcmp(motorModel.dataSet.TypeOfRotor,'EESM')
    Ir = motorModel.FluxMap_dq.Ir;
    Fr = motorModel.FluxMap_dq.Fr;
end

%% Surfaces
figNames{1} = 'CurrentD';
figNames{2} = 'CurrentQ';
figNames{3} = 'Torque';

for ii=1:length(figNames)
    hfig(ii) = figure();
    figSetting();
    hax(ii) = axes('OuterPosition',[0 0 1 1],...
        'XLim',[min(Fd(:)) max(Fd(:))],...
        'YLim',[min(Fq(:)) max(Fq(:))],...
        'PlotBoxAspectRatio',[1 1 0.8]);
    xlabel('$\lambda_d$ (Vs)')
    ylabel('$\lambda_q$ (Vs)')
    view(3)
    switch ii
        case 1
            zlabel('$i_d$ (A)')
            set(gca,'ZLim',[min(Id(:)) max(Id(:))])
        case 2
            zlabel('$i_q$ (A)')
            set(gca,'ZLim',[min(Iq(:)) max(Iq(:))])
        case 3
            zlabel('$T$ (Nm)')
            set(gca,'ZLim',[min(T(:)) max(T(:))])
    end
    set(hfig(ii),'FileName',[pathname resFolder figNames{ii} '.fig'])
    set(hfig(ii),'Name',figNames{ii})
end

if ~strcmp(motorModel.dataSet.TypeOfRotor,'EESM')
    surf(hax(1),Fd,Fq,Id,'FaceColor','interp','EdgeColor','none')
    contour3(hax(1),Fd,Fq,Id,'EdgeColor','k','ShowText','off')
    surf(hax(2),Fd,Fq,Iq,'FaceColor','interp','EdgeColor','none')
    contour3(hax(2),Fd,Fq,Iq,'EdgeColor','k','ShowText','off')
    surf(hax(3),Fd,Fq,T,'FaceColor','interp','EdgeColor','none')
    contour3(hax(3),Fd,Fq,T,'EdgeColor','k','ShowText','off')
else
    surf(hax(1),Fd(:,:,end),Fq(:,:,end),Id(:,:,end),'FaceColor','interp','EdgeColor','none')
    contour3(hax(1),Fd(:,:,end),Fq(:,:,end),Id(:,:,end),'EdgeColor','k','ShowText','off')
    surf(hax(2),Fd(:,:,end),Fq(:,:,end),Iq(:,:,end),'FaceColor','interp','EdgeColor','none')
    contour3(hax(2),Fd(:,:,end),Fq(:,:,end),Iq(:,:,end),'EdgeColor','k','ShowText','off')
    surf(hax(3),Fd(:,:,end),Fq(:,:,end),T(:,:,end),'FaceColor','interp','EdgeColor','none')
    contour3(hax(3),Fd(:,:,end),Fq(:,:,end),T(:,:,end),'EdgeColor','k','ShowText','off')
    % Multi Ir plot
    size_Ir = size(Ir(1,1,:));
    Ir_plot_index = round(linspace(1, size_Ir(end), 10)); % 10 è il # di plot tra da vedere
    Ir_plot_index(1) = 1;
    Ir_plot_index(end) = max(size_Ir);
    Ir_plot_index = unique(Ir_plot_index); %Eventuali indici ripetuti vengono rimossi
    Ir_values = Ir(1,1,:);
    Ir_plot_values = Ir_values(Ir_plot_index);
    hfig(4) = surfplot_slider(Fd(:,:,end),Fq(:,:,end),T(:,:,Ir_plot_index),Ir_plot_values,[pathname resFolder 'Torque_Ir.fig'],[{'$\lambda_d$ (Vs)'}, {'$\lambda_q$ (Vs)'}, {'T (Nm)'}]);
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
