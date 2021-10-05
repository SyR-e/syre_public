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
resFolder = [motName '_results\MMM results\' 'Inverse Model dq - ' int2str(motorModel.data.tempPM) 'deg\'];

% Load data
Id = motorModel.FluxMapInv_dq.Id;
Iq = motorModel.FluxMapInv_dq.Iq;
Fd = motorModel.FluxMapInv_dq.Fd;
Fq = motorModel.FluxMapInv_dq.Fq;
T  = motorModel.FluxMapInv_dq.T;

%% Surfaces
figNames{1} = 'CurrentD';
figNames{2} = 'CurrentQ';
figNames{3} = 'Torque';

for ii=1:length(figNames)
    hfig(ii) = figure();
    figSetting();
    hax(ii) = axes('OuterPosition',[0 0 1 1],...
        'XLim',[min(min(Fd)) max(max(Fd))],...
        'YLim',[min(min(Fq)) max(max(Fq))],...
        'PlotBoxAspectRatio',[1 1 0.8]);
    xlabel('$\lambda_d$ [$Vs$]')
    ylabel('$\lambda_q$ [$Vs$]')
    view(3)
    switch ii
        case 1
            zlabel('$i_d$ [$A$]')
            set(gca,'ZLim',[min(min(Id)) max(max(Id))])
        case 2
            zlabel('$i_q$ [$A$]')
            set(gca,'ZLim',[min(min(Iq)) max(max(Iq))])
        case 3
            zlabel('$T$ [$Nm$]')
            set(gca,'ZLim',[min(min(T)) max(max(T))])
    end
    set(hfig(ii),'FileName',[pathname resFolder figNames{ii} '.fig'])
end

surf(hax(1),Fd,Fq,Id,'FaceColor','interp','EdgeColor','interp')
surf(hax(2),Fd,Fq,Iq,'FaceColor','interp','EdgeColor','interp')
surf(hax(3),Fd,Fq,T,'FaceColor','interp','EdgeColor','interp')

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
