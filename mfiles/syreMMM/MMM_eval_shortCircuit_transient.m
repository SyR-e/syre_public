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

function [TSCout,resFolderOut] = MMM_eval_shortCircuit_transient(motorModel,saveFlag)

if nargin()==1
    saveFlag = [];
end

% load data
fdfq     = motorModel.FluxMap_dq;
ironLoss = motorModel.IronPMLossMap_dq;
ACloss   = motorModel.acLossFactor;
pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;

p           = motorModel.data.p;
n3phase     = motorModel.data.n3phase;
kActiveSets = 1;

iAmp   = motorModel.WaveformSetup.CurrAmpl;
gamma  = motorModel.WaveformSetup.CurrAngle;
n      = motorModel.WaveformSetup.EvalSpeed;
nCycle = motorModel.WaveformSetup.nCycle;

nPoints = 1000;

if strcmp(motorModel.WaveformSetup.IronLossFlag,'Yes')
    ironLossFactor = motorModel.WaveformSetup.IronLossFactor;
    if strcmp(motorModel.WaveformSetup.PMLossFlag,'Yes')
        PMLossFactor = motorModel.WaveformSetup.PMLossFactor;
    else
        PMLossFactor = 0;
    end
    [Pfe,~,~,~,~,Ppm] = calcIronLoss(ironLoss,fdfq,n*p/60);
    Pfe = Pfe*ironLossFactor+Ppm*PMLossFactor;
else
    ironLossFactor = 0;
    PMLossFactor = 0;
    Pfe = zeros(size(fdfq.Id));
end

if strcmp(motorModel.WaveformSetup.ACLossFlag,'Yes')
    Rs0    = motorModel.data.Rs;
    tempCu = motorModel.data.tempCu;
    l      = motorModel.data.l;
    lend   = motorModel.data.lend;
    [Rs,~] = calcRsTempFreq(Rs0,tempCu,l,lend,ACloss,'LUT',tempCu,n*p/60);
else
    Rs0    = motorModel.data.Rs;
    Rs     = Rs0;
    tempCu = motorModel.data.tempCu;
end

if ~isempty(motorModel.DemagnetizationLimit)
    Idemag = interp1(motorModel.DemagnetizationLimit.tempPM,motorModel.DemagnetizationLimit.Idemag,motorModel.data.tempPM);
else
    Idemag = NaN;
end


if isfield(motorModel.WaveformSetup,'flagPlot')
    flagPlot = motorModel.WaveformSetup.flagPlot;
else
    if nCycle<3
        flagPlot = 1;
    else
        flagPlot = 0;
    end
end

if n==0
    error('Speed cannot be 0!!!')
end

% mirror the flux maps, if needed
if strcmp(motorModel.data.axisType,'SR')
    if min(fdfq.Id,[],'all')>=0
        fdfq.Id = [-fliplr(fdfq.Id(:,2:end)) fdfq.Id];
        fdfq.Iq = [fdfq.Iq(:,2:end) fdfq.Iq];
        fdfq.Fd = [-fliplr(fdfq.Fd(:,2:end)) fdfq.Fd];
        fdfq.Fq = [fliplr(fdfq.Fq(:,2:end)) fdfq.Fq];
        fdfq.T  = [-fliplr(fdfq.T(:,2:end)) fdfq.T];
        Pfe     = [+fliplr(Pfe(:,2:end)) Pfe];
        % fdfq.We = [+fliplr(fdfq.We(:,2:end)) fdfq.We];
        % fdfq.Wc = [+fliplr(fdfq.Wc(:,2:end)) fdfq.Wc];
    end
    if strcmp(motorModel.data.motorType,'SR')
        if min(fdfq.Iq,[],'all')>=0
            fdfq.Id = [fdfq.Id(2:end,:);fdfq.Id];
            fdfq.Iq = [-flipud(fdfq.Iq(2:end,:));fdfq.Iq];
            fdfq.Fd = [flipud(fdfq.Fd(2:end,:));fdfq.Fd];
            fdfq.Fq = [-flipud(fdfq.Fq(2:end,:));fdfq.Fq];
            fdfq.T  = [-flipud(fdfq.T(2:end,:));fdfq.T];
            Pfe     = [+flipud(Pfe(2:end,:));Pfe];
            % fdfq.We = [+flipud(fdfq.We(2:end,:)) fdfq.We];
            % fdfq.Wc = [+flipud(fdfq.Wc(2:end,:)) fdfq.Wc];
        end
    end
else
    if min(fdfq.Iq,[],'all')>=0
        fdfq.Id = [fdfq.Id(2:end,:);fdfq.Id];
        fdfq.Iq = [-flipud(fdfq.Iq(2:end,:));fdfq.Iq];
        fdfq.Fd = [flipud(fdfq.Fd(2:end,:));fdfq.Fd];
        fdfq.Fq = [-flipud(fdfq.Fq(2:end,:));fdfq.Fq];
        fdfq.T  = [-flipud(fdfq.T(2:end,:));fdfq.T];
        Pfe     = [+flipud(Pfe(2:end,:));Pfe];
        % fdfq.We = [+flipud(fdfq.We(2:end,:)) fdfq.We];
        % fdfq.Wc = [+flipud(fdfq.Wc(2:end,:)) fdfq.Wc];
    end
end

% prepare maps
IdqM = fdfq.Id+j*fdfq.Iq;
Fdq  = fdfq.Fd+j*fdfq.Fq;
Tem  = fdfq.T;

% prepare iron loss maps
Vind  = 1j*2*pi*(n*p/60).*(fdfq.Fd+j*fdfq.Fq);
IdqFe = 2/3*Pfe./conj(Vind);
IdqFe(Pfe==0) = 0;
%Rfe = 3/2*abs(n*pi/30*p*(fdfq.Fd+j*fdfq.Fq)).^2./Pfe;
IdqFe = smoothdata2(real(IdqFe))+j*smoothdata2(imag(IdqFe));

Tm = Tem-Pfe./(n*pi/30);



%% initial values
w = n*pi/30*p;
t1 = 2*pi/w; % electrical period
time = linspace(0,nCycle*t1,nPoints*nCycle+1);

idVect      = nan(size(time));
iqVect      = nan(size(time));
fdVect      = nan(size(time));
fqVect      = nan(size(time));
TemVect     = nan(size(time));
TemVectExt  = nan(size(time));
TmVect      = nan(size(time));
TmVectExt   = nan(size(time));
intVect     = ones(size(time));
idmVect     = nan(size(time));
iqmVect     = nan(size(time));
idFeVect    = nan(size(time));
iqFeVect    = nan(size(time));
pCuDCVect   = nan(size(time));
pCuACVect   = nan(size(time));
pCuVect     = nan(size(time));
pFeVect     = nan(size(time));
% weVect      = nan(size(time));
% wcVect      = nan(size(time));

id0 = iAmp.*cosd(gamma);
iq0 = iAmp.*sind(gamma);
fd0 = interp2(fdfq.Id,fdfq.Iq,fdfq.Fd,id0,iq0);
fq0 = interp2(fdfq.Id,fdfq.Iq,fdfq.Fq,id0,iq0);

idmVect(1)    = iAmp*cosd(gamma);
iqmVect(1)    = iAmp*sind(gamma);
fdVect(1)     = interp2(fdfq.Id,fdfq.Iq,fdfq.Fd,id0,iq0);
fqVect(1)     = interp2(fdfq.Id,fdfq.Iq,fdfq.Fq,id0,iq0);
TemVect(1)    = interp2(fdfq.Id,fdfq.Iq,Tem,id0,iq0);
TemVectExt(1) = TemVect(1);
TmVect(1)     = interp2(fdfq.Id,fdfq.Iq,Tm,id0,iq0);
TmVectExt(1)  = TmVect(1);
idFeVect(1)   = interp2(fdfq.Id,fdfq.Iq,real(IdqFe),id0,iq0);
iqFeVect(1)   = interp2(fdfq.Id,fdfq.Iq,imag(IdqFe),id0,iq0);
idVect(1)     = idmVect(1)+idFeVect(1);
iqVect(1)     = iqmVect(1)+iqFeVect(1);
pCuDCVect(1)  = 3/2*n3phase*kActiveSets*Rs0*abs(idVect(1)+j*iqVect(1)).^2;
pCuVect(1)    = 3/2*n3phase*kActiveSets*Rs*abs(idVect(1)+j*iqVect(1)).^2;
pCuACVect(1)  = pCuVect(1)-pCuDCVect(1);
pFeVect(1)    = interp2(fdfq.Id,fdfq.Iq,Pfe,id0,iq0);
% weVect(1)     = interp2(fdfq.Id,fdfq.Iq,fdfq.We,id0,iq0);
% wcVect(1)     = interp2(fdfq.Id,fdfq.Iq,fdfq.Wc,id0,iq0);

%% steady-state solution
debug.nMax   = n;
debug.tempCu = tempCu;
[SSSCout,~] = MMM_eval_shortCircuit_steadyState(motorModel,0,debug);

%% interpolant functions (inverse maps, magnetic, 4Q)
data.Fd  = fdfq.Fd(:);
data.Fq  = fdfq.Fq(:);
data.Id  = fdfq.Id(:);
data.Iq  = fdfq.Iq(:);
data.Tem = fdfq.T(:);
data.Tm  = Tm(:);
data.IdFe = real(IdqFe(:));
data.IqFe = real(IdqFe(:));

filt = ones(1,length(data.Id));
filt(isnan(data.Fd)) = 0;
filt(isnan(data.Fq)) = 0;

data.Fd(filt==0)  = [];
data.Fq(filt==0)  = [];
data.Id(filt==0)  = [];
data.Iq(filt==0)  = [];
data.Tem(filt==0) = [];
data.Tm(filt==0)  = [];

fInt.Id  = scatteredInterpolant(data.Fd,data.Fq,data.Id,'linear','none');
fInt.Iq  = scatteredInterpolant(data.Fd,data.Fq,data.Iq,'linear','none');
fInt.Tem = scatteredInterpolant(data.Fd,data.Fq,data.Tem,'linear','none');
fInt.Tm  = scatteredInterpolant(data.Fd,data.Fq,data.Tm,'linear','none');


%% prepare figures
iStr = num2str(abs(id0+j*iq0),3);
iStr = strrep(iStr,'.','A');
if ~contains(iStr,'A')
    iStr=[iStr 'A'];
end
gammaStr=num2str(atan2(iq0,id0)*180/pi,4);
gammaStr=strrep(gammaStr,'.','d');
if ~contains(gammaStr,'d')
    gammaStr=[gammaStr 'd'];
end
resFolder = [motName '_results\MMM results\' 'TransientSC - ' iStr '_' gammaStr '_' int2str(n) 'rpm' '_' int2str(motorModel.data.tempPM) 'degPM_' int2str(motorModel.data.tempCu) 'degCu\'];

resFolderOut = [pathname resFolder];

hfig(1) = figure();
figSetting();
xlim([min(time) max(time)]);
hax(1) = gca;
xlabel('$t$ (s)')
ylabel('(A)')
hleg(1) = legend(hax(1),'show','Location','southoutside');
set(gcf,'FileName',[pathname resFolder 'CurrentVStime.fig'])
plot(hax(1),[time(1) time(end)],real(SSSCout.Im)*[1 1],'--b','DisplayName','$i_{d,m}$ - steady-state')
plot(hax(1),[time(1) time(end)],imag(SSSCout.Im)*[1 1],'--r','DisplayName','$i_{q,m}$ - steady-state')
plot(hax(1),[time(1) time(end)],abs(SSSCout.Im)*[1 1],'--k','DisplayName','$|i_{m}|$ - steady-state')
hp1(1) = plot(hax(1),time,idVect.*intVect,':b','DisplayName','$i_{d,m}$','HandleVisibility','off');
hp1(2) = plot(hax(1),time,iqVect.*intVect,':r','DisplayName','$i_{q,m}$','HandleVisibility','off');
hp1(3) = plot(hax(1),time,abs(idVect+j*iqVect).*intVect,':k','DisplayName','$|i_{m}|$','HandleVisibility','off');
hp1(4) = plot(hax(1),time,idVect,'-b','DisplayName','$i_{d,m}$');
hp1(5) = plot(hax(1),time,iqVect,'-r','DisplayName','$i_{q,m}$');
hp1(6) = plot(hax(1),time,abs(idVect+j*iqVect),'-k','DisplayName','$|i_{m}|$');
if ~isnan(Idemag)
    plot(hax(1),[time(1) time(end)],Idemag*[1 1],'-','Color',[0.5 0 0],'LineWidth',2,'DisplayName','$I_{demag}$');
end


hfig(2) = figure();
figSetting();
set(gca,...
    'DataAspectRatio',[1 1 1],...
    'GridAlpha',1,...
    'Layer','top');
xlabel('$i_{d,m}$ (A)')
ylabel('$i_{q,m}$ (A)')
title('Flux linkage during short circuit and current trajectory')
set(gcf,'FileName',[pathname resFolder 'idiqTrajectory.fig'])
[c,h] = contourf(fdfq.Id,fdfq.Iq,abs(fdfq.Fd+j*fdfq.Fq),'DisplayName','$|\lambda_{dq}|$');
clabel(c,h);
tmp = get(gca,'CLim');
set(gca,'CLim',tmp,'CLimMode','manual')
[c,h] = contour(fdfq.Id,fdfq.Iq,abs(fdfq.Fd+j*fdfq.Fq),abs(fd0+j*fq0)*[1 1],'-k','LineWidth',1.5,'DisplayName','$|\lambda_{ini}|$');
[c,h] = contour(fdfq.Id,fdfq.Iq,abs(fdfq.Id+j*fdfq.Iq),abs(id0+j*iq0)*[1 1],'--k','LineWidth',1.5,'DisplayName','$|i_{m,ini}|$');
hax(2) = gca;
hp2(1) = plot(hax(2),idVect,iqVect,'-r','DisplayName','Transient');
plot(hax(2),real(SSSCout.Im),imag(SSSCout.Im),'r*','DisplayName','Steady-state')
hleg(2) = legend(hax(2),'show','Location','northeast');

hfig(3) = figure();
figSetting();
xlim([min(time) max(time)]);
hax(3) = gca;
xlabel('$t$ (s)')
ylabel('$T$ (Nm)')
hleg(3) = legend(hax(3),'show','Location','southeast');
set(gcf,'FileName',[pathname resFolder 'TorqueVStime.fig'])
plot(hax(3),[time(1) time(end)],SSSCout.Tem*[1 1],'--b','DisplayName','Steady-state')
hp3(1) = plot(hax(3),time,TemVect,'-b','DisplayName','$T_{em}$ - Interpolated');
hp3(2) = plot(hax(3),time,TemVectExt,':b','DisplayName','$T_{em}$ - Extrapolated');
hp3(3) = plot(hax(3),time,TmVect,'-r','DisplayName','$T_{m}$ - Interpolated');
hp3(4) = plot(hax(3),time,TmVectExt,':r','DisplayName','$T_{m}$ - Extrapolated');

hfig(4) = figure();
figSetting();
xlim([min(time) max(time)]);
hax(4) = gca;
xlabel('$t$ (s)')
ylabel('$\lambda$ (Vs)')
hleg(4) = legend(hax(4),'show','Location','southeast');
set(gcf,'FileName',[pathname resFolder 'FluxVStime.fig'])
plot(hax(4),[time(1) time(end)],real(SSSCout.Fdq)*[1 1],'--b','DisplayName','$\lambda_d$ - steady-state')
plot(hax(4),[time(1) time(end)],imag(SSSCout.Fdq)*[1 1],'--r','DisplayName','$\lambda_q$ - steady-state')
plot(hax(4),[time(1) time(end)],abs(SSSCout.Fdq)*[1 1],'--k','DisplayName','$|\lambda_{dq}|$ - steady-state')
hp4(1) = plot(hax(4),time,fdVect.*intVect,':b','DisplayName','$\lambda_d$','HandleVisibility','off');
hp4(2) = plot(hax(4),time,fqVect.*intVect,':r','DisplayName','$\lambda_q$','HandleVisibility','off');
hp4(3) = plot(hax(4),time,abs(fdVect+j*fqVect).*intVect,':k','DisplayName','$|\lambda_{dq}|$','HandleVisibility','off');
hp4(4) = plot(hax(4),time,fdVect,'-b','DisplayName','$\lambda_d$');
hp4(5) = plot(hax(4),time,fqVect,'-r','DisplayName','$\lambda_q$');
hp4(6) = plot(hax(4),time,abs(fdVect+j*fqVect),'-k','DisplayName','$|\lambda_{dq}|$');

FdPlot = linspace(min(data.Fd,[],'all'),max(data.Fd,[],'all'),256);
FqPlot = linspace(min(data.Fq,[],'all'),max(data.Fq,[],'all'),256);
[FdPlot,FqPlot] = meshgrid(FdPlot,FqPlot);
IdPlot = fInt.Id(FdPlot,FqPlot);
IqPlot = fInt.Iq(FdPlot,FqPlot);

hfig(5) = figure();
figSetting();
set(gca,...
    'DataAspectRatio',[1 1 1],...
    'GridAlpha',1,...
    'Layer','top');
xlabel('$\lambda_d$ (Vs)')
ylabel('$\lambda_q$ (Vs)')
title('Current during short circuit and flux linkage trajectory')
set(gcf,'FileName',[pathname resFolder 'fdfqTrajectory.fig'])
[c,h] = contourf(FdPlot,FqPlot,abs(IdPlot+j*IqPlot),'DisplayName','$|i_{m}|$');
clabel(c,h);
[c,h] = contour(FdPlot,FqPlot,abs(IdPlot+j*IqPlot),abs(id0+j*iq0)*[1 1],'-k','LineWidth',1.5,'DisplayName','$|i_{ini}|$');
tmp = get(gca,'CLim');
set(gca,'CLim',tmp,'CLimMode','manual')
[c,h] = contour(FdPlot,FqPlot,abs(FdPlot+j*FqPlot),abs(fdVect(1)+j*fqVect(1))*[1 1],'--k','LineWidth',1.5,'DisplayName','$|\lambda_{ini}|$');
hax(5) = gca;
hp5(1) = plot(hax(5),fdVect,fqVect,'-r','DisplayName','Transient');
plot(hax(5),real(SSSCout.Fdq),imag(SSSCout.Fdq),'r*','DisplayName','Steady-state')
hleg(5) = legend(hax(5),'show','Location','northeast');

plot(hax(2),idVect(1),iqVect(1),'ro','HandleVisibility','off');
plot(hax(5),fdVect(1),fqVect(1),'ro','HandleVisibility','off');

hfig(6) = figure();
figSetting();
xlim([min(time) max(time)]);
hax(6) = gca;
xlabel('$t$ (s)')
ylabel('$P_{loss}$ (W)')
hleg(6) = legend(hax(6),'show','Location','northeast');
set(gcf,'FileName',[pathname resFolder 'LossVStime.fig'])
tmp = [pCuDCVect;pCuACVect;pFeVect];
hp6 = area(hax(6),time,tmp');
set(hp6(1),'DisplayName','$P_{Cu,DC}$');
set(hp6(2),'DisplayName','$P_{Cu,AC}$');
set(hp6(3),'DisplayName','$P_{Fe}$');

hfig(7) = figure();
figSetting();
xlim([min(time) max(time)]);
hax(7) = gca;
xlabel('$t$ (s)')
ylabel('$|I_{dq}|$ (A)')
hleg(7) = legend(hax(7),'show','Location','northeast');
set(gcf,'FileName',[pathname resFolder 'CurrentComponentsVStime.fig'])
hp7(1) = plot(hax(7),time,abs(idmVect+j*iqmVect),'-','DisplayName','$i_m$');
hp7(2) = plot(hax(7),time,abs(idFeVect+j*iqFeVect),'-','DisplayName','$i_{Fe}$');
hp7(3) = plot(hax(7),time,abs(idVect+j*iqVect),'-','DisplayName','$i_{ph}$');

for ii=1:length(hleg)
    switch ii
        case {3,4,5}
            set(hleg(ii),'NumColumns',2,'Location','southoutside');
        case {1,2}
            set(hleg(ii),'NumColumns',2,'Location','southoutside');
    end
end

for ii=1:length(hfig)
    tmp = get(hfig(ii),'FileName');
    [~,name,~] = fileparts(tmp);
    set(hfig(ii),'Name',name);
end


%% solution
fInt.Id   = scatteredInterpolant(data.Fd,data.Fq,data.Id,'linear','linear');
fInt.Iq   = scatteredInterpolant(data.Fd,data.Fq,data.Iq,'linear','linear');
fInt.Tem  = scatteredInterpolant(data.Fd,data.Fq,data.Tem,'linear','linear');
fInt.Tm   = scatteredInterpolant(data.Fd,data.Fq,data.Tm,'linear','linear');
fInt.IdFe = scatteredInterpolant(data.Fd,data.Fq,real(IdqFe(:)),'linear','linear');
fInt.IqFe = scatteredInterpolant(data.Fd,data.Fq,imag(IdqFe(:)),'linear','linear');
dt = time(2);

disp('Transient short-circuit computation...')
fprintf(' %06.2f%%',0)

for tt=2:length(time)
    dFd = (+w*fqVect(tt-1)-n3phase*Rs*idVect(tt-1))*dt;
    dFq = (-w*fdVect(tt-1)-n3phase*Rs*iqVect(tt-1))*dt;
    
    fdVect(tt)    = fdVect(tt-1)+dFd;
    fqVect(tt)    = fqVect(tt-1)+dFq;
    
    idmVect(tt)   = fInt.Id(fdVect(tt),fqVect(tt));
    iqmVect(tt)   = fInt.Iq(fdVect(tt),fqVect(tt));
    TemVect(tt)   = interp2(fdfq.Id,fdfq.Iq,Tem,idmVect(tt),iqmVect(tt));
    TmVect(tt)    = interp2(fdfq.Id,fdfq.Iq,Tm,idmVect(tt),iqmVect(tt));
    idFeVect(tt)  = fInt.IdFe(fdVect(tt),fqVect(tt));
    iqFeVect(tt)  = fInt.IqFe(fdVect(tt),fqVect(tt));
    % idFeVect(tt)  = interp2(fdfq.Id,fdfq.Iq,real(IdqFe),idmVect(tt),iqmVect(tt));
    % iqFeVect(tt)  = interp2(fdfq.Id,fdfq.Iq,imag(IdqFe),idmVect(tt),iqmVect(tt));

    if sign(idmVect(tt))==sign(idFeVect(tt))
        idFeVect(tt) = -idFeVect(tt);
    end
    if sign(iqmVect(tt))==sign(iqFeVect(tt))
        iqFeVect(tt) = -iqFeVect(tt);
    end
    if abs(idFeVect(tt))>abs(idmVect(tt))
        idFeVect(tt) = sign(idFeVect(tt))*abs(idmVect(tt));
    end
    if abs(iqFeVect(tt))>abs(iqmVect(tt))
        iqFeVect(tt) = sign(iqFeVect(tt))*abs(iqmVect(tt));
    end

    idVect(tt)     = idmVect(tt)+idFeVect(tt);
    iqVect(tt)     = iqmVect(tt)+iqFeVect(tt);
    TemVectExt(tt) = fInt.Tem(fdVect(tt),fqVect(tt));
    TmVectExt(tt)  = fInt.Tm(fdVect(tt),fqVect(tt));
    intVect = ~isnan(TemVect);
    intVect = double(intVect);
    intVect(intVect==0) = NaN;

    pCuDCVect(tt) = 3/2*n3phase*kActiveSets*Rs0*abs(idVect(tt)+j*iqVect(tt))^2;
    pCuVect(tt)   = 3/2*n3phase*kActiveSets*Rs*abs(idVect(tt)+j*iqVect(tt))^2;
    pCuACVect(tt)  = pCuVect(tt)-pCuDCVect(tt);
    pFeVect(tt)    = interp2(fdfq.Id,fdfq.Iq,Pfe,idmVect(tt),iqmVect(tt));
    % weVect(tt)     = interp2(fdfq.Id,fdfq.Iq,fdfq.We,idmVect(tt),iqmVect(tt));
    % wcVect(tt)     = interp2(fdfq.Id,fdfq.Iq,fdfq.Wc,idmVect(tt),iqmVect(tt));

    fprintf('\b\b\b\b\b\b\b\b')
    fprintf(' %06.2f%%',tt/length(time)*100)

    if flagPlot
        set(hp1(1),'YData',idVect);
        set(hp1(2),'YData',iqVect);
        set(hp1(3),'YData',abs(idVect+j*iqVect));
        set(hp1(4),'YData',idVect.*intVect);
        set(hp1(5),'YData',iqVect.*intVect);
        set(hp1(6),'YData',abs(idVect+j*iqVect).*intVect);
        
        set(hp2(1),'XData',idVect,'YData',iqVect);

        set(hp3(1),'YData',TemVect);
        set(hp3(2),'YData',TemVectExt);
        set(hp3(3),'YData',TmVect);
        set(hp3(4),'YData',TmVectExt);
        
        set(hp4(1),'YData',fdVect);
        set(hp4(2),'YData',fqVect);
        set(hp4(3),'YData',abs(fdVect+j*fqVect));
        set(hp4(4),'YData',fdVect.*intVect);
        set(hp4(5),'YData',fqVect.*intVect);
        set(hp4(6),'YData',abs(fdVect+j*fqVect).*intVect);
        
        set(hp5(1),'XData',fdVect,'YData',fqVect);

        set(hp6(1),'YData',pCuDCVect);
        set(hp6(2),'YData',pCuACVect);
        set(hp6(3),'YData',pFeVect);

        set(hp7(1),'YData',abs(idmVect+j*iqmVect));
        set(hp7(2),'YData',abs(idFeVect+j*iqFeVect));
        set(hp7(3),'YData',abs(idVect+j*iqVect));
        
        drawnow()
    end
end

disp(' ')
disp('Transient short-circuit computed!')

if ~flagPlot
    set(hp1(1),'YData',idVect);
    set(hp1(2),'YData',iqVect);
    set(hp1(3),'YData',abs(idVect+j*iqVect));
    set(hp1(4),'YData',idVect.*intVect);
    set(hp1(5),'YData',iqVect.*intVect);
    set(hp1(6),'YData',abs(idVect+j*iqVect).*intVect);

    set(hp2(1),'XData',idVect,'YData',iqVect);

    set(hp3(1),'YData',TemVect);
    set(hp3(2),'YData',TemVectExt);
    set(hp3(3),'YData',TmVect);
    set(hp3(4),'YData',TmVectExt);

    set(hp4(1),'YData',fdVect);
    set(hp4(2),'YData',fqVect);
    set(hp4(3),'YData',abs(fdVect+j*fqVect));
    set(hp4(4),'YData',fdVect.*intVect);
    set(hp4(5),'YData',fqVect.*intVect);
    set(hp4(6),'YData',abs(fdVect+j*fqVect).*intVect);

    set(hp5(1),'XData',fdVect,'YData',fqVect);

    set(hp6(1),'YData',pCuDCVect);
    set(hp6(2),'YData',pCuACVect);
    set(hp6(3),'YData',pFeVect);

    set(hp7(1),'YData',abs(idmVect+j*iqmVect));
    set(hp7(2),'YData',abs(idFeVect+j*iqFeVect));
    set(hp7(3),'YData',abs(idVect+j*iqVect));

    drawnow()
end

%% output data
TSCout.time    = time;
TSCout.iph     = idVect+j*iqVect;
TSCout.im      = idmVect+j*iqmVect;
TSCout.ife     = idFeVect+j*iqFeVect;
TSCout.Tem     = TemVectExt;
TSCout.Tm      = TmVectExt;
TSCout.intFilt = intVect;
TSCout.fdq     = fdVect+j*fqVect;
TSCout.n       = n;
TSCout.SSSC    = SSSCout;
TSCout.Rs      = Rs0*ones(size(time));
TSCout.RsDC    = Rs*ones(size(time));
TSCout.pCu     = pCuVect;
TSCout.pCuDC   = pCuDCVect;
TSCout.pCuAC   = pCuACVect;
TSCout.pFe     = pFeVect;
% ASCout.we      = weVect;
% ASCout.wc      = wcVect;

% ASCout.fmap.IdqM  = IdqM;
% ASCout.fmap.IdqFe = IdqFe;
% ASCout.fmap.Idq   = IdqM+IdqFe;
% ASCout.fmap.Fdq   = Fdq;
% ASCout.fmap.Vdq   = Vdq;
% ASCout.fmap.Pfe   = Pfe;
% ASCout.fmap.Tem   = Tem;
% ASCout.fmap.Tm    = Tm;


%% save data and figures

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
    
    save([pathname resFolder 'TransientShortCircuitResults.mat'],'motorModel','TSCout');

    for ii=1:length(hfig)
        savePrintFigure(hfig(ii));
    end
end
