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

function MMM_plot_FFTdqtMap(motorModel)

dqtMap    = motorModel.FluxMap_dqt;
hVect     = motorModel.dqtElab.harmonic;
pathname  = motorModel.data.pathname;
motName   = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'Torque Ripple FFT analysis - ' int2str(motorModel.data.tempPM) 'deg\'];
Id        = dqtMap.data.Id;
Iq        = dqtMap.data.Iq;

% computation
a = fft(dqtMap.data.T,[],3);
aL = size(a,3);
FFT.cont = 1*abs(a(:,:,1))/aL;
FFT.harm = 2*abs(a(:,:,2:floor(aL/2)))/aL;

% figures
for ii=1:length(hVect)
    hfig(ii) = figure();
    figSetting();
    hax(ii) = axes('OuterPosition',[0 0 1 1],...
        'DataAspectRatio',[1 1 1],...
        'Layer','top',...
        'XLim',[min(Id,[],'all') max(Id,[],'all')],...
        'YLim',[min(Iq,[],'all') max(Iq,[],'all')]);
    xlabel('$i_d$ [$A$]')
    ylabel('$i_q$ [$A$]')
    title(['$\Delta T_{' int2str(hVect(ii)) '}$']);
    set(hfig(ii),'FileName',[pathname resFolder 'Torque Ripple Harmonic ' int2str(hVect(ii)) '.fig']);
    [c,h] = contourf(dqtMap.data.Id(:,:,1),dqtMap.data.Iq(:,:,1),FFT.harm(:,:,hVect(ii)),'DisplayName',['$\Delta T_{' int2str(hVect(ii)) '}$']);
    clabel(c,h)
end

% save figure
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




