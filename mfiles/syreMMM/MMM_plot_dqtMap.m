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

function MMM_plot_dqtMap(motorModel)

% load data
dqtMap = motorModel.FluxMap_dqt;

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'dqtMap Model - ' int2str(motorModel.data.tempPM) 'deg\'];

if ~exist([pathname resFolder],'dir')
    mkdir([pathname resFolder]);
end

IdLim = [min(dqtMap.data.Id,[],'all') max(dqtMap.data.Id,[],'all')];
IqLim = [min(dqtMap.data.Iq,[],'all') max(dqtMap.data.Iq,[],'all')];
FdLim = [min(dqtMap.data.Fd,[],'all') max(dqtMap.data.Fd,[],'all')];
FqLim = [min(dqtMap.data.Fq,[],'all') max(dqtMap.data.Fq,[],'all')];
TLim  = [min(dqtMap.data.T,[],'all') max(dqtMap.data.T,[],'all')];

% flux d
figure()
figSetting()
view(3)
set(gca,...
    'XLim',IdLim,...
    'YLim',IqLim,...
    'ZLim',FdLim,...
    'CLim',FdLim,...
    'PlotBoxAspectRatio',[1 1 0.8])
xlabel('$i_d$ [$A$]')
ylabel('$i_q$ [$A$]')
zlabel('$\lambda_d$ [$Vs$]')
plotData.x  = dqtMap.data.Id;
plotData.y  = dqtMap.data.Iq;
plotData.z  = dqtMap.data.Fd;
plotData.th = dqtMap.th;
plot_dqtMap(gca,plotData,[pathname resFolder 'FluxD.gif'])

% flux q
figure()
figSetting()
view(3)
set(gca,...
    'XLim',IdLim,...
    'YLim',IqLim,...
    'ZLim',FqLim,...
    'CLim',FqLim,...
    'PlotBoxAspectRatio',[1 1 0.8])
xlabel('$i_d$ [$A$]')
ylabel('$i_q$ [$A$]')
zlabel('$\lambda_q$ [$Vs$]')
plotData.x  = dqtMap.data.Id;
plotData.y  = dqtMap.data.Iq;
plotData.z  = dqtMap.data.Fq;
plotData.th = dqtMap.th;
plot_dqtMap(gca,plotData,[pathname resFolder 'FluxQ.gif'])

% torque
figure()
figSetting()
view(3)
set(gca,...
    'XLim',IdLim,...
    'YLim',IqLim,...
    'ZLim',TLim,...
    'CLim',TLim,...
    'PlotBoxAspectRatio',[1 1 0.8])
xlabel('$i_d$ [$A$]')
ylabel('$i_q$ [$A$]')
zlabel('$T$ [$Nm$]')
plotData.x  = dqtMap.data.Id;
plotData.y  = dqtMap.data.Iq;
plotData.z  = dqtMap.data.T;
plotData.th = dqtMap.th;
plot_dqtMap(gca,plotData,[pathname resFolder 'Torque.gif'])




