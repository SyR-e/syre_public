% Copyright 2022
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

function MMM_plot_AppInductanceMap(motorModel)

% load data
Id = motorModel.AppInductanceMap_dq.Id;
Iq = motorModel.AppInductanceMap_dq.Iq;
Ld = motorModel.AppInductanceMap_dq.Ld;
Lq = motorModel.AppInductanceMap_dq.Lq;
Fm = motorModel.AppInductanceMap_dq.Fm;

axisType = motorModel.data.axisType;

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'Apparent Inductance Maps - ' int2str(motorModel.data.tempPM) 'deg\'];

%% Surfaces
figNames{1} = 'InducD';
figNames{2} = 'InducQ';
figNames{3} = 'FluxM';
figNames{4} = 'Anisotropy';

switch axisType
    case 'SR'
        csi = Ld./Lq;
        csiName = '$L_d/L_q$';
    case 'PM'
        csi = Lq./Ld;
        csiName = '$L_q/L_d$';
end

for ii=1:length(figNames)
    hfig(ii) = figure();
    figSetting();
    hax(ii) = axes('OuterPosition',[0 0 1 1],...
        'XLim',[min(Id,[],'all') max(Id,[],'all')],...
        'YLim',[min(Iq,[],'all') max(Iq,[],'all')],...
        'PlotBoxAspectRatio',[1 1 0.8]);
    xlabel('$i_d$ [A]')
    ylabel('$i_q$ [A]')
    view(3)
    set(hfig(ii),'FileName',[pathname resFolder figNames{ii} '.fig'])
    set(hfig(ii),'Name',figNames{ii});
    switch ii
        case 1
            zlabel('$L_{d}$ [H]')
            set(gca,'ZLim',[min(Ld,[],'all') max(Ld,[],'all')]);
        case 2
            zlabel('$L_{q}$ [H]')
            set(gca,'ZLim',[min(Lq,[],'all') max(Lq,[],'all')]);
        case 3
            zlabel('$\lambda_m$ [Vs]')
            set(gca,'ZLim',[min(Fm,[],'all') max(Fm,[],'all')]);
        case 4
            zlabel(csiName)
            set(gca,'ZLim',[min(csi,[],'all') max(csi,[],'all')]);
    end
end

surf(hax(1),Id,Iq,Ld,'FaceColor','interp','EdgeColor','interp')
surf(hax(2),Id,Iq,Lq,'FaceColor','interp','EdgeColor','interp')
surf(hax(3),Id,Iq,Fm,'FaceColor','interp','EdgeColor','interp')
surf(hax(4),Id,Iq,csi,'FaceColor','interp','EdgeColor','interp')


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


