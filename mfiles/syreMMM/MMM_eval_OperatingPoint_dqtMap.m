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

dqtMap    = motorModel.FluxMap_dqt;
iAmp      = motorModel.WaveformSetup.CurrAmpl;
gamma     = motorModel.WaveformSetup.CurrAngle;
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
    resFolder=[motorName '_results\MMM results\' 'T_eval_' iStr '_' gammaStr ' - ' int2str(motorModel.data.tempPM) 'deg\'];
    mkdir(pathname,resFolder);
    SOL.th = dqtMap.th;
    SOL.id = id(ii)*ones(size(SOL.th));
    SOL.iq = iq(ii)*ones(size(SOL.th));
    % SOL.fd = dqtMap.fInt.Fd(SOL.id,SOL.iq,SOL.th);
    % SOL.fq = dqtMap.fInt.Fq(SOL.id,SOL.iq,SOL.th);
    % SOL.T  = dqtMap.fInt.T(SOL.id,SOL.iq,SOL.th);
    SOL.fd = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fd,SOL.id,SOL.iq,SOL.th,'spline');
    SOL.fq = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fq,SOL.id,SOL.iq,SOL.th,'spline');
    SOL.T  = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.T,SOL.id,SOL.iq,SOL.th,'spline');

    if isfield(dqtMap.data,'Fa')
        % SOL.fa = dqtMap.fInt.Fa(SOL.id,SOL.iq,SOL.th);
        % SOL.fb = dqtMap.fInt.Fb(SOL.id,SOL.iq,SOL.th);
        % SOL.fc = dqtMap.fInt.Fc(SOL.id,SOL.iq,SOL.th);
        SOL.fa = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fa,SOL.id,SOL.iq,SOL.th,'spline');
        SOL.fb = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fb,SOL.id,SOL.iq,SOL.th,'spline');
        SOL.fc = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fc,SOL.id,SOL.iq,SOL.th,'spline');
    else
        SOL.fa = nan(size(SOL.th));
        SOL.fb = nan(size(SOL.th));
        SOL.fc = nan(size(SOL.th));
    end
    if isfield(dqtMap.data,'Ia')
        % SOL.ia = dqtMap.fInt.Ia(SOL.id,SOL.iq,SOL.th);
        % SOL.ib = dqtMap.fInt.Ib(SOL.id,SOL.iq,SOL.th);
        % SOL.ic = dqtMap.fInt.Ic(SOL.id,SOL.iq,SOL.th);
        SOL.ia = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Ia,SOL.id,SOL.iq,SOL.th,'spline');
        SOL.ib = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Ib,SOL.id,SOL.iq,SOL.th,'spline');
        SOL.ic = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Ic,SOL.id,SOL.iq,SOL.th,'spline');
    else
        SOL.ia = nan(size(SOL.th));
        SOL.ib = nan(size(SOL.th));
        SOL.ic = nan(size(SOL.th));
    end
    
    if isfield(dqtMap,'sets')
        n3phase = length(dqtMap.sets);
        for ss=1:n3phase
            SOL.ia(ss,:) = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Ia,SOL.id,SOL.iq,SOL.th,'spline');
            SOL.ib(ss,:) = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Ib,SOL.id,SOL.iq,SOL.th,'spline');
            SOL.ic(ss,:) = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Ic,SOL.id,SOL.iq,SOL.th,'spline');
            SOL.fa(ss,:) = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Fa,SOL.id,SOL.iq,SOL.th,'spline');
            SOL.fb(ss,:) = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Fb,SOL.id,SOL.iq,SOL.th,'spline');
            SOL.fc(ss,:) = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Fc,SOL.id,SOL.iq,SOL.th,'spline');
            
            SOL.sets(ss).id = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Id,SOL.id,SOL.iq,SOL.th,'spline');
            SOL.sets(ss).iq = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Iq,SOL.id,SOL.iq,SOL.th,'spline');
            SOL.sets(ss).fd = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Fd,SOL.id,SOL.iq,SOL.th,'spline');
            SOL.sets(ss).fq = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Fq,SOL.id,SOL.iq,SOL.th,'spline');

            SOL.sets(ss).ia = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Ia,SOL.id,SOL.iq,SOL.th,'spline');
            SOL.sets(ss).ib = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Ib,SOL.id,SOL.iq,SOL.th,'spline');
            SOL.sets(ss).ic = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Ic,SOL.id,SOL.iq,SOL.th,'spline');
            SOL.sets(ss).fa = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Fa,SOL.id,SOL.iq,SOL.th,'spline');
            SOL.sets(ss).fb = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Fb,SOL.id,SOL.iq,SOL.th,'spline');
            SOL.sets(ss).fc = interpn(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.sets(ss).Fc,SOL.id,SOL.iq,SOL.th,'spline');
        end
    end
    
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
    
%     if ~isempty(motorModel.IronPMLossMap_dq)
%         fdfq = motorModel.FluxMap_dq;
%         ironLoss = motorModel.IronPMLossMap_dq;
%         FreqElet=motorModel.WaveformSetup.EvalSpeed/60*motorModel.data.p;
%         [~,Pfes_h,Pfes_c,Pfer_h,Pfer_c,Ppm] = calcIronLoss(ironLoss,fdfq,FreqElet);
%         out.Pfes_h = interp2(fdfq.Id,fdfq.Iq,Pfes_h,id(ii),iq(ii));
%         out.Pfes_c = interp2(fdfq.Id,fdfq.Iq,Pfes_c,id(ii),iq(ii));
%         out.Pfer_h = interp2(fdfq.Id,fdfq.Iq,Pfer_h,id(ii),iq(ii));
%         out.Pfer_c = interp2(fdfq.Id,fdfq.Iq,Pfer_c,id(ii),iq(ii));
%         out.Ppm    = interp2(fdfq.Id,fdfq.Iq,Ppm,id(ii),iq(ii));
%         out.velDim = motorModel.WaveformSetup.EvalSpeed;
%     end
    
%     save([pathname resFolder motorName '_T_eval_' iStr '_' gammaStr '.mat'],'out');
    save([pathname resFolder 'out.mat'],'out','motorModel');
    
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
    
    figure()
    figSetting()
    subplot(2,1,1)
    set(gca,'XLim',[0 360],'XTick',0:60:360);
    xlabel('$\theta$ [elt deg]')
    ylabel('$\lambda_{abc}$ [Vs]')
    title(['Phase flux linkages'])
    plot(SOL.th,SOL.fa);
    plot(SOL.th,SOL.fb);
    plot(SOL.th,SOL.fc);
    subplot(2,1,2)
    set(gca,'XLim',[0 360],'XTick',0:60:360);
    xlabel('$\theta$ [elt deg]')
    ylabel('$i_{abc}$ [A]')
    title(['Phase currents'])
    plot(SOL.th,SOL.ia);
    plot(SOL.th,SOL.ib);
    plot(SOL.th,SOL.ic);
    
    saveas(gcf,[pathname resFolder 'singt_' iStr '_' gammaStr '_plot_Phase.fig'])
    
    output{ii}=out;
end

if length(id)>1
    resFolder=[motorName '_results\MMM results\senseOut - ' int2str(motorModel.data.tempPM) 'deg - ' datestr(now,30) '\'];
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





