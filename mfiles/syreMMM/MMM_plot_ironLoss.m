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

function MMM_plot_ironLoss(motorModel)

fdfq     = motorModel.FluxMap_dq;
ironLoss = motorModel.IronPMLossMap_dq;

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'Iron Loss Model - ' int2str(motorModel.data.tempPM) 'deg\'];

n0 = ironLoss.n0;
f0 = ironLoss.f0;

figNames{1} = 'StatorHysteresis';
figNames{2} = 'StatorEddyCurrents';
figNames{3} = 'RotorHysteresis';
figNames{4} = 'RotorEddyCurrents';
figNames{5} = 'PermanentMagnets';
figNames{6} = 'IronLoss';

nPU = [0.5 1 2];

colors{1} = [1.0 0.0 0.0];
colors{2} = [0.0 0.0 1.0];
colors{3} = [0.0 0.8 0.0];

for ii=1:length(figNames)
    hfig(ii) = figure();
    figSetting()
    hax(ii) = axes('OuterPosition',[0 0 1 1],...
        'XLim',[min(fdfq.Id,[],'all') max(fdfq.Id,[],'all')],...
        'YLim',[min(fdfq.Iq,[],'all') max(fdfq.Iq,[],'all')],...
        'PlotBoxAspectRatio',[1 1 0.8]);
    view(3)
    xlabel('$i_d$ (A)')
    ylabel('$i_q$ (A)')
    switch ii
        case 1
            zlabel('$P_{Fe,s,h}$ (W)')
            title('Stator Hysteresis')
        case 2
            zlabel('$P_{Fe,s,c}$ (W)')
            title('Stator Eddy-currents')
        case 3
            zlabel('$P_{Fe,r,h}$ (W)')
            title('Rotor Hysteresis')
        case 4
            zlabel('$P_{Fe,r,c}$ (W)')
            title('Rotor Eddy-currents')
        case 5
            zlabel('$P_{PM}$ (W)')
            title('Permanent Magnets')
        case 6
            zlabel('$P_{Fe}$ (W)')
            title('Total Iron Loss')
    end
    set(hfig(ii),'FileName',[pathname resFolder figNames{ii} '.fig'])
    set(hfig(ii),'Name',figNames{ii})
end


for ii=1:length(nPU)
    [Pfe,Pfesh,Pfesc,Pferh,Pferc,Ppm] = calcIronLoss(ironLoss,fdfq,f0*nPU(ii));
    surf(hax(1),fdfq.Id,fdfq.Iq,Pfesh,'FaceColor',colors{ii},'EdgeColor',0.5*colors{ii},'DisplayName',['$n=' int2str(n0*nPU(ii)) '$ rpm'])
    surf(hax(2),fdfq.Id,fdfq.Iq,Pfesc,'FaceColor',colors{ii},'EdgeColor',0.5*colors{ii},'DisplayName',['$n=' int2str(n0*nPU(ii)) '$ rpm'])
    surf(hax(3),fdfq.Id,fdfq.Iq,Pferh,'FaceColor',colors{ii},'EdgeColor',0.5*colors{ii},'DisplayName',['$n=' int2str(n0*nPU(ii)) '$ rpm'])
    surf(hax(4),fdfq.Id,fdfq.Iq,Pferc,'FaceColor',colors{ii},'EdgeColor',0.5*colors{ii},'DisplayName',['$n=' int2str(n0*nPU(ii)) '$ rpm'])
    surf(hax(5),fdfq.Id,fdfq.Iq,Ppm,'FaceColor',colors{ii},'EdgeColor',0.5*colors{ii},'DisplayName',['$n=' int2str(n0*nPU(ii)) '$ rpm'])
    surf(hax(6),fdfq.Id,fdfq.Iq,Pfe,'FaceColor',colors{ii},'EdgeColor',0.5*colors{ii},'DisplayName',['$n=' int2str(n0*nPU(ii)) '$ rpm'])
end

for ii=1:length(hax)
    hleg(ii) = legend(hax(ii),'show','Location','northeast');
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







