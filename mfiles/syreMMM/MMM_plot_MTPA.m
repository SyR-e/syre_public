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

function MMM_plot_MTPA(motorModel)

MTPA = motorModel.AOA.MTPA;
MTPV = motorModel.AOA.MTPV;
fdfq = motorModel.fdfq;
pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'AOA - ' int2str(motorModel.data.tempPM) 'deg\'];

%% Figures
% torque and torque ripple VS current
hfig(1) = figure();
figSetting()
xlabel('$I$ [$A$]')
ylabel('$T$ [$Nm$]')
title('Torque vs peak current along MTPA')
set(hfig(1),'FileName',[pathname resFolder 'torqueVScurrent.fig'])
plot(abs(MTPA.id+j*MTPA.iq),MTPA.T,'-b','DisplayName','$T$')
plot(abs(MTPA.id+j*MTPA.iq),MTPA.T+MTPA.dTpp/2,'-r','LineWidth',1)
plot(abs(MTPA.id+j*MTPA.iq),MTPA.T-MTPA.dTpp/2,'-r','LineWidth',1)

hfig(2) = figure();
figSetting()
xlabel('$T$ [$Nm$]')
ylabel('$\Delta T_{pp}$ [$Nm$]')
title('Torque ripple vs torque along the MTPA')
set(hfig(2),'FileName',[pathname resFolder 'rippleVStorque.fig'])
plot(MTPA.T,MTPA.dTpp,'-b')

hfig(3) = figure();
figSetting()
xlabel('$I$ [$A$]')
ylabel('$k_T$ [$Nm/A$]')
title('Torque constant vs peak current along MTPA')
set(hfig(3),'FileName',[pathname resFolder 'ktVScurrent.fig'])
plot(abs(MTPA.id(2:end)+j*MTPA.iq(2:end)),MTPA.T(2:end)./abs(MTPA.id(2:end)+j*MTPA.iq(2:end)),'-b')

hfig(4) = figure();
figSetting();
xlabel('$I$ [$A$]')
ylabel('$\gamma$ [$^\circ$]')
title('Current angle along the MTPA')
set(hfig(4),'FileName',[pathname resFolder 'gammaVScurrent.fig'])
plot(abs(MTPA.id+j*MTPA.iq),angle(MTPA.id+j*MTPA.iq)*180/pi,'-b')

hfig(5) = figure();
figSetting()
xlabel('$T_{set}$ [$Nm$]')
ylabel('$i_{dq}$ [$A$]')
title('MTPA LUT')
set(hfig(5),'FileName',[pathname resFolder 'tablesMTPA.fig'])
plot(MTPA.T,MTPA.id,'-b','DisplayName','$i_d$ set-point')
plot(MTPA.T,MTPA.iq,'-r','DisplayName','$i_q$ set-point')
% plot(MTPA.T,MTPA.id,'kx','DisplayName','$i_d$ LUT')
% plot(MTPA.T,MTPA.iq,'kx','DisplayName','$i_q$ LUT')
hleg = legend('show','Location','northwest');
set(hleg,'NumColumns',2)

hfig(6) = figure();
figSetting()
xlabel('$\lambda_{MTPV}$ [$Vs$]')
ylabel('$T_{MTPV}$ [$Nm$]')
title('MTPV LUT')
set(hfig(6),'FileName',[pathname resFolder 'tablesMTPV.fig'])
plot(abs(MTPV.fd+j*MTPV.fq),MTPV.T,'-b','DisplayName','curve')
% plot(abs(MTPV.fd+j*MTPV.fq),MTPV.T,'kx','DisplayName','LUT')
hleg = legend('show','Location','northwest');

hfig(7) = figure();
figSetting();
xlabel('$T$ [$Nm$]')
ylabel('$\delta$ [$^\circ$]')
title('Flux angle versus torque')
set(hfig(7),'FileName',[pathname resFolder 'deltaVStorque.fig'])
plot(MTPA.T,angle(MTPA.fd+j*MTPA.fq)*180/pi,'-b','DisplayName','MTPA')
plot(MTPV.T,angle(MTPV.fd+j*MTPV.fq)*180/pi,'-r','DisplayName','MTPV')
hleg = legend('show','Location','southeast');

hfig(8) = figure();
figSetting();
xlabel('$T$ [$Nm$]')
ylabel('$\lambda$ [$Vs$]')
title('Flux amplitude versus torque')
set(hfig(8),'FileName',[pathname resFolder 'fluxVStorque.fig'])
plot(MTPA.T,abs(MTPA.fd+j*MTPA.fq),'-b','DisplayName','MTPA')
plot(MTPV.T,abs(MTPV.fd+j*MTPV.fq),'-r','DisplayName','MTPV')
hleg = legend('show','Location','southeast');

hfig(9) = figure();
figSetting();
set(gca,...
    'XLim',[min(fdfq.Id,[],'all') max(fdfq.Id,[],'all')],...
    'YLim',[min(fdfq.Iq,[],'all') max(fdfq.Iq,[],'all')],...
    'DataAspectRatio',[1 1 1]);
xlabel('$i_d$ [$A$]')
ylabel('$i_q$ [$A$]')
iLevels = get(gca,'XTick');
set(gca,'XTick',iLevels,'YTick',iLevels);
set(hfig(9),'FileName',[pathname resFolder 'MTPAMTPV(id,iq).fig'])
[c,h] = contour(fdfq.Id,fdfq.Iq,fdfq.T,'-','LineColor',[0.0 0.8 0.0],'LineWidth',1,'DisplayName','$T$ [$Nm$]');
clabel(c,h)
[c,h] = contour(fdfq.Id,fdfq.Iq,abs(fdfq.Id+j*fdfq.Iq),iLevels,'-','LineColor',[0.0 0.0 0.0],'LineWidth',0.5,'DisplayName','$I$ [$A$]');
clabel(c,h)
plot(MTPA.id,MTPA.iq,'-r','DisplayName','MTPA')
plot(MTPV.id,MTPV.iq,'-b','DisplayName','MTPV')
hleg = legend('show','Location','northeast');


hfig(10) = figure();
figSetting();
set(gca,...
    'XLim',[min(fdfq.Id,[],'all') max(fdfq.Id,[],'all')],...
    'YLim',[min(fdfq.Iq,[],'all') max(fdfq.Iq,[],'all')],...
    'DataAspectRatio',[1 1 1]);
if ~sum(isnan(fdfq.dTpp))
    set(gca,'Clim',[min(fdfq.dTpp,[],'all') max(fdfq.dTpp,[],'all')]);
else
    set(gca,'Clim',[0 1]);
end
xlabel('$i_d$ [$A$]')
ylabel('$i_q$ [$A$]')
iLevels = get(gca,'XTick');
set(gca,'XTick',iLevels,'YTick',iLevels);
set(hfig(10),'FileName',[pathname resFolder 'TorqueRippleMap.fig'])
if ~sum(isnan(fdfq.dTpp))
    [c,h] = contourf(fdfq.Id,fdfq.Iq,fdfq.dTpp,'LineWidth',1,'DisplayName','$\Delta T_{pp}$ [$Nm$]');
    clabel(c,h)
else
    [c,h] = contourf(fdfq.Id,fdfq.Iq,zeros(size(fdfq.dTpp)),'LineWidth',1,'DisplayName','$\Delta T_{pp}$ [$Nm$]');
    clabel(c,h)
end

[c,h] = contour(fdfq.Id,fdfq.Iq,fdfq.T,'-','LineColor',[0.0 0.0 0.0],'LineWidth',1,'DisplayName','$T$ [$Nm$]');
clabel(c,h)
[c,h] = contour(fdfq.Id,fdfq.Iq,abs(fdfq.Id+j*fdfq.Iq),iLevels,'-','LineColor',[1.0 0.0 0.0],'LineWidth',0.5,'DisplayName','$I$ [$A$]');
clabel(c,h)
plot(MTPA.id,MTPA.iq,'-r','DisplayName','MTPA')
plot(MTPV.id,MTPV.iq,'-r','DisplayName','MTPV')
hleg = legend('show','Location','northeast');


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