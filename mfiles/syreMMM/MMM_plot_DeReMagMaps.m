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

function MMM_plot_DeReMagMaps(motorModel,saveFlag)

DFVC = motorModel.VFMdata.Control.DFVC;
FOC = motorModel.VFMdata.Control.FOC;
fdfq = motorModel.FluxMap_dq;
pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = checkPathSyntax([motName '_results\MMM results\' 'DeReMagMapsControl - ' int2str(motorModel.data.tempPM) 'deg\']); %% chiedere come cambairlo

%% Figures

%% FOC
    % Demag Id vs Ms and Torque
    hfig(1) = figure();
    figSetting()
    xlabel('$MS$')
    ylabel('$T$ (Nm)')
    zlabel('$Id$ (Nm)')
    title('Demag Id')
    set(hfig(1),'FileName',[pathname resFolder 'RemagId.fig'])
    surf(FOC.Demag.MS, FOC.Demag.T, FOC.Demag.id)
    
    % Demag Iq vs Ms and Torque
    hfig(2) = figure();
    figSetting()
    xlabel('$MS$')
    ylabel('$T$ (Nm)')
    zlabel('$Iq$ (Nm)')
    title('Demag Iq')
    set(hfig(2),'FileName',[pathname resFolder 'DemagIq.fig'])
    surf(FOC.Demag.MS, FOC.Demag.T, FOC.Demag.iq)
    
    % 1-D LUT Torque Limit FOC demag 
    hfig(3) = figure();
    figSetting()
    xlabel('$MS$',Interpreter='latex')
    ylabel('$T_{max}$ (Nm)')
    title('Max Torque in FOC during Demag')
    set(hfig(3),'FileName',[pathname resFolder 'Demag1DTorqueLimitFOC.fig'])
    plot(FOC.Demag.Limit.MSlim, FOC.Demag.Limit.Tlim, '-b')

    % Remag Id vs Ms and Torque
    hfig(4) = figure();
    figSetting()
    xlabel('$MS$')
    ylabel('$T$ (Nm)')
    zlabel('$Id$ (Nm)')
    title('Remag Id')
    set(hfig(4),'FileName',[pathname resFolder 'RemagId.fig'])
    surf(FOC.Remag.MS, FOC.Remag.T, FOC.Remag.id)
    
    % Remag Iq vs Ms and Torque
    hfig(5) = figure();
    figSetting()
    xlabel('$MS$')
    ylabel('$T$ (Nm)')
    zlabel('$Iq$ (Nm)')
    title('Remag Iq')
    set(hfig(5),'FileName',[pathname resFolder 'RemagIq.fig'])
    surf(FOC.Remag.MS, FOC.Remag.T, FOC.Remag.iq)
    
    % 1-D LUT Torque Limit FOC remag 
    hfig(6) = figure();
    figSetting()
    xlabel('$MS$')
    ylabel('$T_{max}$ (Nm)')
    title('Max Torque in FOC during Remag')
    set(hfig(6),'FileName',[pathname resFolder 'Remag1DTorqueLimitFOC.fig'])
    plot(FOC.Remag.Limit.MSlim, FOC.Remag.Limit.Tlim, '-b')
    
%% DFVC    
    % Demag Flux vs Ms and Torque 
    hfig(7) = figure();
    figSetting()
    xlabel('$MS$')
    ylabel('$T$ (Nm)')
    zlabel('$F$ (Wb)')
    title('Demag Flux')
    set(hfig(7),'FileName',[pathname resFolder 'DemagFlx.fig'])
    surf(DFVC.Demag.MS, DFVC.Demag.T, DFVC.Demag.F)
    
    % 1-D LUT Torque Limit DFVC demag 
    hfig(8) = figure();
    figSetting()
    xlabel('$MS$ ',Interpreter='latex')
    ylabel('$T_{max}$ (Nm)')
    title('Max Torque in DFVC during Demag')
    set(hfig(8),'FileName',[pathname resFolder 'Demag1DTorqueLimitDFVC.fig'])
    plot(DFVC.Demag.Limit.MSlim, DFVC.Demag.Limit.Tlim, '-b')
    
    % Remag Flux vs Ms and Torque 
    hfig(9) = figure();
    figSetting()
    xlabel('$MS$')
    ylabel('$T$ (Nm)')
    zlabel('$F$ (Wb)')
    title('Remag Flux')
    set(hfig(9),'FileName',[pathname resFolder 'DemagFlx.fig'])
    surf(DFVC.Demag.MS, DFVC.Demag.T, DFVC.Demag.F)
    
    % 1-D LUT Torque Limit DFVC remag 
    hfig(10) = figure();
    figSetting()
    xlabel('$MS$')
    ylabel('$T_{max}$ (Nm)')
    title('Max Torque in DFVC during Remag')
    set(hfig(10),'FileName',[pathname resFolder 'Remag1DTorqueLimitDFVC.fig'])
    plot(DFVC.Remag.Limit.MSlim, DFVC.Remag.Limit.Tlim, '-b')

for ii=1:length(hfig)
    tmp = get(hfig(ii),'FileName');
    [~,name,~] = fileparts(tmp);
    set(hfig(ii),'Name',name);
end


%% Save figures
if nargin()==1
    answer = 'No';
    answer = questdlg('Save figures?','Save','Yes','No',answer);
else
    if saveFlag
        answer = 'Yes';
    else
        answer = 'No';
    end
end

if strcmp(answer,'Yes')
    if ~exist([pathname resFolder],'dir')
        mkdir([pathname resFolder]);
    end
    
    for ii=1:length(hfig)
        savePrintFigure(hfig(ii));
    end
end