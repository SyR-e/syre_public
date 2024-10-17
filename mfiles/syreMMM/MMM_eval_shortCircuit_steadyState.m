% Copyright 2024
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

function [SSSCout,resFolderOut] = MMM_eval_shortCircuit_steadyState(motorModel,saveFlag,debug)

figFlag = 1;

if nargin()==1
    saveFlag = [];
    debug    = [];
elseif nargin()==2
    debug = [];
end

if ~isempty(debug)
    %saveFlag = 0;
    figFlag = 0;
end

fdfq         = motorModel.FluxMap_dq;
ironLoss     = motorModel.IronPMLossMap_dq;
AClossFactor = motorModel.acLossFactor;

p           = motorModel.data.p;
Rs0         = motorModel.data.Rs;
tempCu      = motorModel.data.tempCu;
l           = motorModel.data.l;
lend        = motorModel.data.lend;
n3phase     = motorModel.data.n3phase;
kActiveSets = 1;
nMax        = motorModel.WaveformSetup.EvalSpeed;
nPoints     = 501;

if ~isempty(debug)
    nMax    = debug.nMax;
    nPoints = 1;
    tempCu  = debug.tempCu;
end

nVect = logspace(0,log10(nMax),nPoints);
% nVect = linspace(0,nMax,nPoints);

if strcmp(motorModel.WaveformSetup.IronLossFlag,'Yes')
    ironLossFactor = motorModel.WaveformSetup.IronLossFactor;
    if strcmp(motorModel.WaveformSetup.PMLossFlag,'Yes')
        PMLossFactor = motorModel.WaveformSetup.PMLossFactor;
    else
        PMLossFactor = 0;
    end
else
    ironLossFactor = 0;
    PMLossFactor = 0;
end

if strcmp(motorModel.WaveformSetup.ACLossFlag,'Yes')
    AClossFlag = 1;
else
    AClossFlag = 0;
end

Id = fdfq.Id;
Iq = fdfq.Iq;
Fd = fdfq.Fd;
Fq = fdfq.Fq;
T  = fdfq.T;
% We = fdfq.We;

if strcmp(motorModel.data.axisType,'SR')
    if min(Id,[],'all')>=0
        Id = [-fliplr(Id(:,2:end)) Id];
        Iq = [Iq(:,2:end) Iq];
        Fd = [-fliplr(Fd(:,2:end)) Fd];
        Fq = [fliplr(Fq(:,2:end)) Fq];
        T  = [-fliplr(T(:,2:end)) T];
        % We = [+fliplr(We(:,2:end)) We];
    end
else
    if min(Iq,[],'all')>=0
        Id = [Id(2:end,:);Id];
        Iq = [-flipud(Iq(2:end,:));Iq];
        Fd = [flipud(Fd(2:end,:));Fd];
        Fq = [-flipud(Fq(2:end,:));Fq];
        T  = [-flipud(T(2:end,:));T];
        % We = [-flipud(We(2:end,:));We];
    end
end


idVect = Id(1,:);
iqVect = Iq(:,1);

IdqM = Id+j*Iq;
Fdq  = Fd+j*Fq;

w = nVect*pi/30*p;

Iph_SC = zeros(size(nVect));
T_SC   = zeros(size(nVect));
Im_SC  = zeros(size(nVect));
Ife_SC = zeros(size(nVect));
Pjs_SC = zeros(size(nVect));
Pfe_SC = zeros(size(nVect));
Rs_SC  = zeros(size(nVect));
Vph_SC = zeros(size(nVect));

disp('Steady-state short-circuit computation...')
fprintf(' %06.2f%%',0)

for ii=1:length(nVect)
    freq = nVect(ii)/60*p;
    % AC loss factor
    if AClossFlag
        [Rs,~] = calcRsTempFreq(Rs0,tempCu,l,lend,AClossFactor,'LUT',tempCu,freq);
    else
        Rs = Rs0;
    end
    % iron and PM loss
    if ironLossFactor==0
        Pfe = zeros(size(IdqM));
    else
        [Pfe,~,~,~,~,Ppm] = calcIronLoss(ironLoss,fdfq,freq);
        PfeTmp = Pfe*ironLossFactor+Ppm*PMLossFactor;
        Pfe = zeros(size(IdqM));
        if strcmp(motorModel.data.axisType,'PM')
            PfePos = interp2(fdfq.Id,fdfq.Iq,PfeTmp,real(IdqM),imag(IdqM));
            PfeNeg = interp2(fdfq.Id,-fdfq.Iq,PfeTmp,real(IdqM),imag(IdqM));
            Pfe(imag(IdqM)>=0) = PfePos(imag(IdqM)>=0);
            Pfe(imag(IdqM)<0) = PfeNeg(imag(IdqM)<0);
        else
            PfePos = interp2(fdfq.Id,fdfq.Iq,PfeTmp,real(IdqM),imag(IdqM));
            PfeNeg = interp2(-fdfq.Id,fdfq.Iq,PfeTmp,real(IdqM),imag(IdqM));
            Pfe(real(IdqM)>=0) = PfePos(real(IdqM)>=0);
            Pfe(real(IdqM)<0) = PfeNeg(real(IdqM)<0);
        end
    end
    Pfe(Pfe<0) = 0;
    IdqFe=2/3/(n3phase*kActiveSets)*Pfe./conj(j*w(ii).*Fdq);
    IdqFe(Pfe==0)=0;
    IdFe = real(IdqFe);
    IqFe = imag(IdqFe);
    IdFe = smoothdata2(IdFe);
    IqFe = smoothdata2(IqFe);
    IdqFe = IdFe+j*IqFe;
    
    % SSSC point computation
    Idq = IdqM+IdqFe;
    Vdq = n3phase*Rs*Idq+j*w(ii)*Fdq;
    c = contourc(idVect,iqVect,real(Vdq),[0 0]);
    if c(2,1)==size(c,2)-1
        idTmp = c(1,2:end);
        iqTmp = c(2,2:end);
        iiTmp = 1:1:numel(idTmp);
        vqTmp = interp2(real(IdqM),imag(IdqM),imag(Vdq),idTmp,iqTmp);
        index = interp1(vqTmp,iiTmp,0);
        idOK  = interp1(iiTmp,idTmp,index);
        iqOK  = interp1(iiTmp,iqTmp,index);
    else
        idOK = NaN;
        iqOK = NaN;
    end
    % if isempty(c)
    %     c = contourc(idVect,iqVect,real(Vdq),min(abs(real(Vdq)))*[1 1]);
    % end
    % idTmp = c(1,2:end);
    % iqTmp = c(2,2:end);
    % iiTmp = 1:1:numel(idTmp);
    % vqTmp = interp2(real(IdqM),imag(IdqM),imag(Vdq),idTmp,iqTmp);
    % index = interp1(vqTmp,iiTmp,0);
    % if isempty(index)
    %     index = interp1(vqTmp,iiTmp,min(abs(imag(Vdq))));
    % end
    % idOK  = interp1(iiTmp,idTmp,index);
    % iqOK  = interp1(iiTmp,iqTmp,index);

    % save point
    Im_SC(ii)  = idOK+j*iqOK;
    Ife_SC(ii) = interp2(real(IdqM),imag(IdqM),IdqFe,idOK,iqOK);
    Iph_SC(ii) = Im_SC(ii)+Ife_SC(ii);
    T_SC(ii)   = interp2(real(IdqM),imag(IdqM),T,idOK,iqOK);
    Pjs_SC(ii) = 3/2*Rs*n3phase*kActiveSets*abs(Iph_SC(ii))^2;
    Pfe_SC(ii) = interp2(real(IdqM),imag(IdqM),Pfe,idOK,iqOK);
    Rs_SC(ii)  = Rs;
    Vph_SC(ii) = interp2(Id,Iq,Vdq,idOK,iqOK);

    fprintf('\b\b\b\b\b\b\b\b')
    fprintf(' %06.2f%%',ii/length(nVect)*100)
    % disp([int2str((ii))])
end

disp(' ')
disp('Steady-state short-circuit computed!')

%output structure

SSSCout.n     = nVect;
SSSCout.Tem   = T_SC;
SSSCout.Iph   = Iph_SC;
SSSCout.Pjs   = Pjs_SC;
SSSCout.PjsDC = 3/2*Rs0*n3phase*kActiveSets*abs(Iph_SC).^2;
SSSCout.PjsAC = SSSCout.Pjs-SSSCout.PjsDC;
SSSCout.Pfe   = Pfe_SC;
SSSCout.Im    = Im_SC;
SSSCout.Ife   = Ife_SC;
SSSCout.Rs    = Rs_SC;
SSSCout.Tm    = SSSCout.Tem-SSSCout.Pfe./(SSSCout.n*pi/30);
SSSCout.Vph   = Vph_SC;
SSSCout.Fdq   = interp2(real(IdqM),imag(IdqM),Fdq,real(Im_SC),imag(Im_SC));

if ~isempty(debug)
    if length(nVect)==1
        SSSCout.fmaps.IdqM  = IdqM;
        SSSCout.fmaps.IdqFe = IdqFe;
        SSSCout.fmaps.Idq   = Idq;
        SSSCout.fmaps.Fdq   = Fdq;
        SSSCout.fmaps.Pfe   = Pfe;
        SSSCout.fmaps.Vdq   = Vdq;
        SSSCout.fmaps.Rs    = Rs;
        SSSCout.fmaps.RsDC  = Rs0;
    end
end

%% outputs
pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;

if figFlag
    resFolder = [motName '_results\MMM results\' 'SteadyStateSC - ' int2str(motorModel.data.tempPM) 'degPM - ' int2str(tempCu) 'degCu' '\'];
    
    resFolderOut = [pathname resFolder];
    
    % create figures
    for ii=1:8
        hfig(ii) = figure();
        figSetting();
        hax(ii) = axes('OuterPosition',[0 0 1 1]);
        switch ii
            case 1
                xlabel('$n$ (rpm)')
                ylabel('$T_{sc}$ (Nm)')
                set(gca,'XLim',[0 nMax])
                set(gcf,'FileName',[pathname resFolder 'torqueVSspeed.fig'])
                legend('show','Location','southeast');
            case 2
                xlabel('$n$ (rpm)')
                ylabel('$I_{sc}$ (A)')
                set(gca,'XLim',[0 nMax]);
                set(gcf,'FileName',[pathname resFolder 'currentVSspeed.fig'])
                legend('show','Location','northeast');
            case 3
                xlabel('$i_{d,m}$ (A)')
                ylabel('$i_{q,m}$ (A)')
                set(gca,...
                    'XLim',[min(Id,[],'all') max(Id,[],'all')],...
                    'YLim',[min(Iq,[],'all') max(Iq,[],'all')],...
                    'DataAspectRatio',[1 1 1],...
                    'CLim',[min(T,[],'all') max(T,[],'all')]);
                set(gcf,'FileName',[pathname resFolder 'DQcurrents.fig']);
                legend('show','Location','northeast');
            case 4
                view(gca,3)
                xlabel('$n$ (rpm)')
                ylabel('$i_d$ (A)')
                zlabel('$i_q$ (A)')
                set(gca,...
                    'XLim',[0 nMax],...
                    'YLim',max(abs(Id),[],'all')*[-1 1],...
                    'ZLim',max(abs(Iq),[],'all')*[-1 1],...
                    'DataAspectRatio',[nMax 2*max(abs(Id),[],'all') 2*max(abs(Iq),[],'all')]);
                set(gcf,'FileName',[pathname resFolder 'idVSiqVSn.fig']);
            case 5
                xlabel('$n$ (rpm)')
                ylabel('$P_{loss}$ (W)')
                set(gca,'XLim',[0 nMax]);
                set(gcf,'FileName',[pathname resFolder 'lossVSspeed.fig'])
                legend('show','Location','northeast');
            case 6
                xlabel('$n$ (rpm)')
                ylabel('$R_s$ ($\Omega$)')
                set(gca,'XLim',[0 nMax]);
                set(gcf,'FileName',[pathname resFolder 'resistanceVSspeed.fig'])
                legend('show','Location','northeast');
            case 7
                xlabel('$n$ (rpm)')
                ylabel('$I_{sc}$ (A)')
                set(gca,'XLim',[0 nMax]);
                set(gcf,'FileName',[pathname resFolder 'currentComponentsVSspeed.fig'])
                legend('show','Location','northeast');
            case 8
                xlabel('$n$ (rpm)')
                ylabel('$P_{loss}$ (W)')
                set(gca,'XLim',[0 nMax]);
                set(gcf,'FileName',[pathname resFolder 'lossVSspeed.fig'])
                legend('show','Location','northeast');
        end
    end
    
    plot(hax(1),nVect,T_SC,'-','DisplayName','$T_{em}$');
    plot(hax(1),nVect,-SSSCout.Pfe./(nVect*pi/30),'-','DisplayName','$T_{Fe}$');
    plot(hax(1),nVect,SSSCout.Tm,'-','DisplayName','$T_{m}$');
    
    plot(hax(2),nVect,real(Iph_SC),'-b','DisplayName','$i_d$')
    plot(hax(2),nVect,imag(Iph_SC),'-r','DisplayName','$i_q$')
    plot(hax(2),nVect,abs(Iph_SC),'--g','DisplayName','$I_{sc}$')
    
    contourf(hax(3),Id,Iq,T,'DisplayName','$T$ [Nm]','ShowText','on');
    contour(hax(3),Id,Iq,abs(Id+j*Iq),'-r','DisplayName','$I$ (A)');
    plot(hax(3),real(Im_SC),imag(Im_SC),'-r','DisplayName','$I_{sc}$ (A)')
    
    plot3(hax(4),nVect,real(Iph_SC),imag(Iph_SC),'-bo');
    
    plot(hax(5),nVect,SSSCout.PjsDC,'DisplayName','$P_{Cu,DC}$')
    plot(hax(5),nVect,SSSCout.PjsAC,'DisplayName','$P_{Cu,AC}$')
    plot(hax(5),nVect,SSSCout.Pfe,'DisplayName','$P_{Fe}$')
    plot(hax(5),nVect,SSSCout.Pjs+SSSCout.Pfe,'-k','DisplayName','$P_{loss}$')
    
    plot(hax(6),nVect,Rs0*ones(size(nVect)),'-','DisplayName','$R_{s,DC}$')
    plot(hax(6),nVect,abs(SSSCout.Rs),'-','DisplayName','$R_s$')
    
    plot(hax(7),nVect,abs(SSSCout.Iph),'-','DisplayName','$I_{ph}$')
    plot(hax(7),nVect,abs(SSSCout.Im),'-','DisplayName','$I_{m}$')
    plot(hax(7),nVect,abs(SSSCout.Ife),'-','DisplayName','$I_{Fe}$')
    
    tmp = [SSSCout.PjsDC;SSSCout.PjsAC;SSSCout.Pfe];
    hPlot = area(hax(8),nVect,tmp');
    
    set(hPlot(1),'DisplayName','$P_{Js,DC}$');
    set(hPlot(2),'DisplayName','$P_{Js,AC}$');
    set(hPlot(3),'DisplayName','$P_{Fe}$');
    
    for ii=1:length(hfig)
        tmp = get(hfig(ii),'FileName');
        [~,name,~] = fileparts(tmp);
        set(hfig(ii),'Name',name);
    end
else
    hfig = [];
    resFolder = [motName '_results\MMM results\' 'SteadyStateSC - ' int2str(motorModel.data.tempPM) 'degPM - ' int2str(tempCu) 'degCu' ' - debug\'];
    resFolderOut = [pathname resFolder];
end
%% Save figures
if isempty(saveFlag)
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

    save([pathname resFolder 'ShortCircuitResults.mat'],'SSSCout')

    for ii=1:length(hfig)
        savePrintFigure(hfig(ii));
    end
end













