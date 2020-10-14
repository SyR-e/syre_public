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

function MMM_eval_OperatingPoint_dqtMap(motorModel)

dqtMap    = motorModel.dqtMap;
iAmp      = motorModel.dqtElab.CurrAmpl;
gamma     = motorModel.dqtElab.CurrAngle;
motorName = motorModel.data.motorName;
pathname  = motorModel.data.pathname;

id = iAmp.*cosd(gamma);
iq = iAmp.*sind(gamma);

for ii=1:length(id)
    iStr=num2str(abs(id(ii)+j*iq(ii)),3);
    iStr=strrep(iStr,'.','A');
    gammaStr=num2str(atan2(iq(ii),id(ii))*180/pi,4);
    gammaStr=strrep(gammaStr,'.','d');
    if ~contains(gammaStr,'d')
        gammaStr=[gammaStr 'd'];
    end
    resFolder=[motorName '_results\MMM results\' 'T_eval' iStr '_' gammaStr ' - ' int2str(motorModel.data.tempPM) 'deg\'];
    mkdir(pathname,resFolder);
    SOL.th = dqtMap.th;
    SOL.id = id(ii)*ones(size(SOL.th));
    SOL.iq = iq(ii)*ones(size(SOL.th));
    SOL.fd = dqtMap.fInt.Fd(SOL.id,SOL.iq,SOL.th);
    SOL.fq = dqtMap.fInt.Fq(SOL.id,SOL.iq,SOL.th);
    SOL.T  = dqtMap.fInt.T(SOL.id,SOL.iq,SOL.th);
    
    out.id   = id(ii);
    out.iq   = iq(ii);
    out.fd   = mean(SOL.fd);
    out.fq   = mean(SOL.fq);
    out.T    = mean(SOL.T);
    out.dT   = std(SOL.T);
    out.dTpu = out.dT/out.T;
    out.dTpp = max(abs(SOL.T))-min(abs(SOL.T));
    out.IPF  = abs(sin(atan2(out.iq,out.id)-atan2(out.fq,out.fd)));
    out.SOL  = SOL;
    
%     save([pathname resFolder motorName '_T_eval_' iStr '_' gammaStr '.mat'],'out');
    save([pathname resFolder 'out.mat'],'out');
    
    figure()
    figSetting()
    subplot(2,1,1)
    plot(SOL.th,SOL.T)
    title(['$T=' num2str(out.T) '\,Nm$']);
    set(gca,'XLim',[0 360],'XTick',0:60:360)
    xlabel('$\theta_{elt}$ [$^\circ$]')
    ylabel('Nm')
    subplot(2,1,2)
    plot(SOL.th,abs(sin(atan2(SOL.iq,SOL.id)-atan2(SOL.fq,SOL.fd))))
    title(['$IPF=' num2str(out.IPF) '$']);
    set(gca,'XLim',[0 360],'XTick',0:60:360)
    xlabel('$\theta_{elt}$ [$^\circ$]')
    ylabel('IPF')
    
    saveas(gcf,[pathname resFolder 'singt_' iStr '_' gammaStr 'plot_Torque.fig']);
    
    figure()
    figSetting()
    subplot(2,1,1)
    plot(SOL.th,SOL.fd)
    title(['$\lambda_d=' num2str(out.fd) 'Vs$'])
    set(gca,'XLim',[0 360],'XTick',0:60:360)
    xlabel('$\theta_{elt}$ [$^\circ$]')
    ylabel('Vs')
    subplot(2,1,2)
    plot(SOL.th,SOL.fq)
    title(['$\lambda_q=' num2str(out.fq) 'Vs$'])
    set(gca,'XLim',[0 360],'XTick',0:60:360)
    xlabel('$\theta_{elt}$ [$^\circ$]')
    ylabel('Vs')
    
    saveas(gcf,[pathname resFolder 'singt_' iStr '_' gammaStr '_plot_Flux.fig'])
    
    output{ii}=out;
end

if length(id)>1
    resFolder=[motorName '_results\MMM results\' motorName '_singT - ' int2str(motorModel.data.tempPM) 'deg\'];
    mkdir(pathname,resFolder);
    %save([pathname resFolder 'senseResults.mat'],'output')
    
    idPlot = zeros(size(id));
    iqPlot = zeros(size(id));
    fdPlot = zeros(size(id));
    fqPlot = zeros(size(id));
    TPlot  = zeros(size(id));
    PFPlot = zeros(size(id));
    dTPlot = zeros(size(id));
    
    for ii=1:length(idPlot)
        idPlot(ii) = output{ii}.id;
        iqPlot(ii) = output{ii}.iq;
        fdPlot(ii) = output{ii}.fd;
        fqPlot(ii) = output{ii}.fq;
        TPlot(ii)  = output{ii}.T;
        PFPlot(ii) = output{ii}.IPF;
        dTPlot(ii) = output{ii}.dTpp;
    end
    
    figure()
    figSetting()
    subplot(2,1,1)
    set(gca,'Xlim',[1 ii],'XTick',1:1:ii);
    xlabel('Simulation \#')
    ylabel('$i_{dq}$ [A]')
    plot(idPlot,'DisplayName','$i_d$')
    plot(iqPlot,'DisplayName','$i_q$')
    plot(abs(idPlot+j*iqPlot),'DisplayName','$|i|$')
    legend('show')
    subplot(2,1,2)
    set(gca,'Xlim',[1 ii],'XTick',1:1:ii);
    xlabel('Simulation \#')
    ylabel('$\gamma$ [$^\circ$]')
    plot(angle(idPlot+j*iqPlot)*180/pi,'DisplayName','$\gamma$')
    legend('show')
    saveas(gcf,[pathname resFolder 'currentSense.fig'])
    
    figure()
    figSetting()
    subplot(2,1,1)
    set(gca,'Xlim',[1 ii],'XTick',1:1:ii);
    xlabel('Simulation \#')
    ylabel('$\lambda_{dq}$ [Vs]')
    plot(fdPlot,'DisplayName','$\lambda_d$')
    plot(fqPlot,'DisplayName','$\lambda_q$')
    plot(abs(fdPlot+j*fqPlot),'DisplayName','$|\lambda|$')
    legend('show')
    subplot(2,1,2)
    set(gca,'Xlim',[1 ii],'XTick',1:1:ii);
    xlabel('Simulation \#')
    ylabel('$\delta$ [$^\circ$]')
    plot(angle(fdPlot+j*fqPlot)*180/pi,'DisplayName','$\delta$')
    legend('show')
    saveas(gcf,[pathname resFolder 'fluxSense.fig'])
    
    figure()
    figSetting()
    subplot(2,1,1)
    set(gca,'Xlim',[1 ii],'XTick',1:1:ii);
    xlabel('Simulation \#')
    ylabel('$T$ [Nm]')
    plot(TPlot,'-')
    plot(TPlot+dTPlot/2,'-r')
    plot(TPlot-dTPlot/2,'-r')
    subplot(2,1,2)
    set(gca,'Xlim',[1 ii],'XTick',1:1:ii);
    xlabel('Simulation \#')
    ylabel('$\Delta T_{pp}$ [Nm]')
    plot(dTPlot)
    saveas(gcf,[pathname resFolder 'torqueSense.fig'])
    
    
    figure()
    figSetting()
    subplot(2,1,1)
    set(gca,'Xlim',[1 ii],'XTick',1:1:ii);
    xlabel('Simulation \#')
    ylabel('$T$ [Nm]')
    plot(TPlot,'-')
    subplot(2,1,2)
    set(gca,'Xlim',[1 ii],'XTick',1:1:ii);
    xlabel('Simulation \#')
    ylabel('$IPF$')
    plot(PFPlot)
    saveas(gcf,[pathname resFolder 'toripfSense.fig'])
    
    senseOut.id   = idPlot;
    senseOut.iq   = iqPlot;
    senseOut.fd   = fdPlot;
    senseOut.fq   = fqPlot;
    senseOut.T    = TPlot;
    senseOut.PF   = PFPlot;
    senseOut.dTpp = dTPlot;
    save([pathname resFolder 'senseResults.mat'],'senseOut')
    
    
end





