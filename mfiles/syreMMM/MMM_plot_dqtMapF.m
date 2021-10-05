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

function MMM_plot_dqtMapF(motorModel)

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'dqtMapF Model - ' int2str(motorModel.data.tempPM) 'deg\'];

if ~exist([pathname resFolder],'dir')
    mkdir([pathname resFolder]);
end

dqtMapF = motorModel.FluxMapInv_dqt;

IdLim = [min(dqtMapF.dataF.Id,[],'all') max(dqtMapF.dataF.Id,[],'all')];
IqLim = [min(dqtMapF.dataF.Iq,[],'all') max(dqtMapF.dataF.Iq,[],'all')];
FdLim = [min(dqtMapF.dataF.Fd,[],'all') max(dqtMapF.dataF.Fd,[],'all')];
FqLim = [min(dqtMapF.dataF.Fq,[],'all') max(dqtMapF.dataF.Fq,[],'all')];
TLim  = [min(dqtMapF.dataF.T,[],'all') max(dqtMapF.dataF.T,[],'all')];

% current d
figure()
figSetting()
view(3)
set(gca,...
    'XLim',FdLim,...
    'YLim',FqLim,...
    'ZLim',IdLim,...
    'CLim',IdLim,...
    'PlotBoxAspectRatio',[1 1 0.8])
xlabel('$\lambda_d$ [$Vs$]')
ylabel('$\lambda_q$ [$Vs$]')
zlabel('$i_d$ [$A$]')
plotData.x  = dqtMapF.dataF.Fd;
plotData.y  = dqtMapF.dataF.Fq;
plotData.z  = dqtMapF.dataF.Id;
plotData.th = dqtMapF.th;
plot_dqtMap(gca,plotData,[pathname resFolder 'CurrentD.gif'])

% current q
figure()
figSetting()
view(3)
set(gca,...
    'XLim',FdLim,...
    'YLim',FqLim,...
    'ZLim',IqLim,...
    'CLim',IqLim,...
    'PlotBoxAspectRatio',[1 1 0.8])
xlabel('$\lambda_d$ [$Vs$]')
ylabel('$\lambda_q$ [$Vs$]')
zlabel('$i_q$ [$A$]')
plotData.x  = dqtMapF.dataF.Fd;
plotData.y  = dqtMapF.dataF.Fq;
plotData.z  = dqtMapF.dataF.Iq;
plotData.th = dqtMapF.th;
plot_dqtMap(gca,plotData,[pathname resFolder 'CurrentQ.gif'])

% torque
figure()
figSetting()
view(3)
set(gca,...
    'XLim',FdLim,...
    'YLim',FqLim,...
    'ZLim',TLim,...
    'CLim',TLim,...
    'PlotBoxAspectRatio',[1 1 0.8])
xlabel('$\lambda_d$ [$Vs$]')
ylabel('$\lambda_q$ [$Vs$]')
zlabel('$T$ [$Nm$]')
plotData.x  = dqtMapF.dataF.Fd;
plotData.y  = dqtMapF.dataF.Fq;
plotData.z  = dqtMapF.dataF.T;
plotData.th = dqtMapF.th;
plot_dqtMap(gca,plotData,[pathname resFolder 'Torque.gif'])

