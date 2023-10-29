% Copyright 2022
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function plot_ASClimit(motorModel,TwMap,demagLimit)

fdfq = motorModel.FluxMap_dq;

TwMap.F = abs(TwMap.Fd+j*TwMap.Fq);

tempPM = motorModel.data.tempPM;

Imax = interp1(demagLimit.temperature,demagLimit.current,tempPM);

switch motorModel.data.axisType
    case 'SR'
        iTmp = unique(fdfq.Iq);
        fTmp = interp2(fdfq.Id,fdfq.Iq,fdfq.Fq,zeros(size(iTmp)),iTmp);
    case 'PM'
        iTmp = unique(fdfq.Id);
        fTmp = -interp2(fdfq.Id,fdfq.Iq,fdfq.Fd,iTmp,zeros(size(iTmp)));
end

    
Fmax = interp1(iTmp,fTmp,Imax);
if isnan(Fmax)
    error('HWC point not found')
end

F0 = interp2(fdfq.Id,fdfq.Iq,abs(fdfq.Fd+j*fdfq.Fq),0,0);
V0 = F0*TwMap.n*pi/30*motorModel.data.p*sqrt(3);
V0(isnan(TwMap.eff)) = NaN;



figure()
figSetting()
xlabel('$n$ [rpm]')
ylabel('$T$ [Nm]')
contourf(TwMap.n,TwMap.T,TwMap.P/1000,'ShowText','on','DisplayName','$P$ [kW]');
plot(unique(TwMap.n),TwMap.T_top_W,'-k','HandleVisibility','off')
plot(unique(TwMap.n),TwMap.T_bot_W,'-k','HandleVisibility','off')
contour(TwMap.n,TwMap.T,TwMap.F,Fmax*[1 1],'-r','LineWidth',2,'DisplayName','ASC limit (HWC)');
legend('show','Location','northeast')
title(['$\Theta_{PM} = ' int2str(tempPM) '^\circ$C'])

figure();
figSetting();
xlabel('$n$ [rpm]')
ylabel('$T$ [Nm]')
plot(TwMap.n(TwMap.F<Fmax),TwMap.T(TwMap.F<Fmax),'g.','DisplayName','ASC safe')
plot(TwMap.n(TwMap.F>Fmax),TwMap.T(TwMap.F>Fmax),'r.','DisplayName','ASC unsafe')
plot(unique(TwMap.n),TwMap.T_top_W,'-k','HandleVisibility','off')
plot(unique(TwMap.n),TwMap.T_bot_W,'-k','HandleVisibility','off')
legend('show','Location','northeast')
title(['Active Short Circuit safe area - $\Theta_{PM} = ' int2str(tempPM) '^\circ$C'])

figure()
figSetting();
xlabel('$n$ [rpm]')
ylabel('$T$ [Nm]')
contourf(TwMap.n,TwMap.T,V0/motorModel.data.Vdc,'ShowText','on','DisplayName','No-load line voltage [p.u.]')
contourf(TwMap.n,TwMap.T,V0/motorModel.data.Vdc,[1 1],'-r','LineWidth',2,'DisplayName','DC link voltage')
plot(unique(TwMap.n),TwMap.T_top_W,'-k','HandleVisibility','off')
plot(unique(TwMap.n),TwMap.T_bot_W,'-k','HandleVisibility','off')
legend('show','Location','northeast')
title(['No-load line voltage - $\Theta_{PM} = ' int2str(tempPM) '^\circ$C'])

figure();
figSetting();
xlabel('$n$ [rpm]')
ylabel('$T$ [Nm]')
plot(TwMap.n(V0<motorModel.data.Vdc),TwMap.T(V0<motorModel.data.Vdc),'g.','DisplayName','Low back-emf')
plot(TwMap.n(V0>motorModel.data.Vdc),TwMap.T(V0>motorModel.data.Vdc),'r.','DisplayName','High back-emf')
plot(unique(TwMap.n),TwMap.T_top_W,'-k','HandleVisibility','off')
plot(unique(TwMap.n),TwMap.T_bot_W,'-k','HandleVisibility','off')
legend('show','Location','northeast')
title(['No-load voltage safe area - $\Theta_{PM} = ' int2str(tempPM) '^\circ$C'])




