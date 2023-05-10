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

function MMM_plot_InductanceMap(motorModel)

% load data
Id  = motorModel.FluxMap_dq.Id;
Iq  = motorModel.FluxMap_dq.Iq;
Ldd = motorModel.IncInductanceMap_dq.Ldd;
Ldq = motorModel.IncInductanceMap_dq.Ldq;
Lqd = motorModel.IncInductanceMap_dq.Lqd;
Lqq = motorModel.IncInductanceMap_dq.Lqq;

axisType = motorModel.data.axisType;

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'Incremental Inductance Maps - ' int2str(motorModel.data.tempPM) 'deg\'];

%% Surfaces
figNames{1} = 'InducDD';
figNames{2} = 'InducDQ';
figNames{3} = 'InducQD';
figNames{4} = 'InducQQ';
figNames{5} = 'Anisotropy';

switch axisType
    case 'SR'
        csi = Ldd./Lqq;
        csiName = '$l_{dd}/l_{qq}$';
    case 'PM'
        csi = Lqq./Ldd;
        csiName = '$l_{qq}/l_{dd}$';
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
    switch ii
        case 1
            zlabel('$l_{dd}$ [H]')
            set(gca,'ZLim',[min(Ldd,[],'all') max(Ldd,[],'all')]);
        case 2
            zlabel('$l_{dq}$ [H]')
            set(gca,'ZLim',[min(Ldq,[],'all') max(Ldq,[],'all')]);
        case 3
            zlabel('$l_{qd}$ [H]')
            set(gca,'ZLim',[min(Lqd,[],'all') max(Lqd,[],'all')]);
        case 4
            zlabel('$l_{qq}$ [H]')
            set(gca,'ZLim',[min(Lqq,[],'all') max(Lqq,[],'all')]);
        case 5
            zlabel(csiName)
            set(gca,'ZLim',[min(csi,[],'all') max(csi,[],'all')]);
    end
end

surf(hax(1),Id,Iq,Ldd,'FaceColor','interp','EdgeColor','interp')
surf(hax(2),Id,Iq,Ldq,'FaceColor','interp','EdgeColor','interp')
surf(hax(3),Id,Iq,Lqd,'FaceColor','interp','EdgeColor','interp')
surf(hax(4),Id,Iq,Lqq,'FaceColor','interp','EdgeColor','interp')
surf(hax(5),Id,Iq,csi,'FaceColor','interp','EdgeColor','interp')


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


