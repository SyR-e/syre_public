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

function MMM_MaxTwPWM_plot(kfix,motorModel,TwMapSIN,TwMapPWM,resFolder)
% resFolder = [motorModel.data.motorName '_results\MMM results\' 'TwMapPWM_' datestr(now,30) '\'];
pathname = motorModel.data.pathname;

n_grid = TwMapPWM.n_grid;
T_grid = TwMapPWM.T_grid;
nmap   = TwMapSIN.nmap;
Tmap   = TwMapSIN.Tmap;
Tmax   = ceil(1.2*max(max(Tmap)/10))*10;

%% Coefficients
%kfe
ii = 1;
hfig(ii)=figure;
figSetting(10,8,8)
set(hfig(ii),'FileName',[pathname resFolder 'PWM_kFe.fig']);
hax(ii) = axes(...
    'XLim',[0 max(max(n_grid))],...
    'YLim',[0 Tmax]);
xlabel('$n$ [rpm]')
ylabel('$T$ [Nm]')
title('PWM $k_{fe}$')
[c,h] = contourf(hax(ii),nmap,Tmap,kfix.kfe);
clabel(c,h)
colorbar

%kjs
ii = ii+1;
hfig(ii)=figure;
figSetting(10,8,8)
set(hfig(ii),'FileName',[pathname resFolder 'PWM_kjs.fig']);
hax(ii) = axes(...
    'XLim',[0 max(max(n_grid))],...
    'YLim',[0 Tmax]);
xlabel('$n$ [rpm]')
ylabel('$T$ [Nm]')
title('PWM $k_{js}$')
[c,h] = contourf(hax(ii),nmap,Tmap,kfix.kjs,1:0.1:max(max(kfix.kjs)));
clabel(c,h)
colorbar

%kPM
ii = ii+1;
hfig(ii)=figure;
figSetting(10,8,8)
set(hfig(ii),'FileName',[pathname resFolder 'PWM_kPM.fig']);
hax(ii) = axes(...
    'XLim',[0 max(max(n_grid))],...
    'YLim',[0 Tmax]);
xlabel('$n$ [rpm]')
ylabel('$T$ [Nm]')
title('PWM $k_{PM}$')
[c,h] = contourf(hax(ii),nmap,Tmap,kfix.kpm);
clabel(c,h)
colorbar

%% Loss 
ii = ii+1;
hfig(ii)=figure;
figSetting(10,8,8)
set(hfig(ii),'FileName',[pathname resFolder 'PWM_IronLoss.fig']);
hax(ii) = axes(...
    'XLim',[0 max(max(n_grid))],...
    'YLim',[0 Tmax]);
xlabel('$n$ [rpm]')
ylabel('$T$ [Nm]')
title('PWM Iron Loss [W]')
[c,h] = contourf(hax(ii),nmap,Tmap,TwMapPWM.Pfe);
clabel(c,h)
colorbar

ii = ii+1;
hfig(ii)=figure;
figSetting(10,8,8)
set(hfig(ii),'FileName',[pathname resFolder 'PWM_CopperLoss.fig']);
hax(ii) = axes(...
    'XLim',[0 max(max(n_grid))],...
    'YLim',[0 Tmax]);
xlabel('$n$ [rpm]')
ylabel('$T$ [Nm]')
title('PWM Copper Loss [W]')
[c,h] = contourf(hax(ii),nmap,Tmap,TwMapPWM.Pjs);
clabel(c,h)
colorbar

ii = ii+1;
hfig(ii)=figure;
figSetting(10,8,8)
set(hfig(ii),'FileName',[pathname resFolder 'PWM_MagnetLoss .fig']);
hax(ii) = axes(...
    'XLim',[0 max(max(n_grid))],...
    'YLim',[0 Tmax]);
xlabel('$n$ [rpm]')
ylabel('$T$ [Nm]')
title('PWM Magnet Loss [W]')
[c,h] = contourf(hax(ii),nmap,Tmap,TwMapPWM.Ppm);
clabel(c,h)
colorbar

%% Effy 
%effy sin
ii = ii+1;
hfig(ii)=figure;
figSetting(10,8,8)
set(hfig(ii),'FileName',[pathname resFolder 'SIN_efficiency.fig']);
hax(ii) = axes(...
    'XLim',[0 max(max(n_grid))],...
    'YLim',[0 Tmax]);
xlabel('$n$ [rpm]')
ylabel('$T$ [Nm]')
title('SIN efficiency')
[c,h] = contourf(hax(ii),nmap,Tmap,TwMapSIN.eff,[64:2:86 87:1:100]/100);
clabel(c,h)
colorbar

%effy pwm
ii = ii+1;
hfig(ii)=figure;
figSetting(10,8,8)
set(hfig(ii),'FileName',[pathname resFolder 'PWM_efficiency.fig']);
hax(ii) = axes(...
    'XLim',[0 max(max(n_grid))],...
    'YLim',[0 Tmax]);
xlabel('$n$ [rpm]')
ylabel('$T$ [Nm]')
title('PWM efficiency')
[c,h] = contourf(hax(ii),nmap,Tmap,TwMapPWM.eff,[64:2:86 87:1:100]/100);
clabel(c,h)
colorbar
plot(hax(ii),n_grid,T_grid,'or')

%% Save figures
for ii=1:length(hfig)
    savePrintFigure(hfig(ii));
end