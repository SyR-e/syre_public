% Copyright 2021
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

function [SCresults] = MMM_eval_shortCircuitTransient(motorModel)

% load data
fdfq   = motorModel.FluxMap_dq;
Rs     = motorModel.data.Rs;
p      = motorModel.data.p;
iAmp   = motorModel.WaveformSetup.CurrAmpl;
gamma  = motorModel.WaveformSetup.CurrAngle;
n      = motorModel.WaveformSetup.EvalSpeed;
nCycle = motorModel.WaveformSetup.nCycle;

nPoints = 1000;
% flagPlot = 1;

if nCycle<3
    flagPlot = 1;
else
    flagPlot = 0;
end

if n==0
    error('Speed cannot be 0!!!')
end


pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;


id0 = iAmp.*cosd(gamma);
iq0 = iAmp.*sind(gamma);

w = n*pi/30*p;

t1 = 2*pi/w; % electrical period

% mirror the flux maps, if needed
if strcmp(motorModel.data.axisType,'SR')
    if min(fdfq.Id,[],'all')>=0
        fdfq.Id = [-fliplr(fdfq.Id(:,2:end)) fdfq.Id];
        fdfq.Iq = [fdfq.Iq(:,2:end) fdfq.Iq];
        fdfq.Fd = [-fliplr(fdfq.Fd(:,2:end)) fdfq.Fd];
        fdfq.Fq = [fliplr(fdfq.Fq(:,2:end)) fdfq.Fq];
        fdfq.T  = [-fliplr(fdfq.T(:,2:end)) fdfq.T];
    end
    if strcmp(motorModel.data.motorType,'SR')
        if min(fdfq.Iq,[],'all')>=0
            fdfq.Id = [fdfq.Id(2:end,:);fdfq.Id];
            fdfq.Iq = [-flipud(fdfq.Iq(2:end,:));fdfq.Iq];
            fdfq.Fd = [flipud(fdfq.Fd(2:end,:));fdfq.Fd];
            fdfq.Fq = [-flipud(fdfq.Fq(2:end,:));fdfq.Fq];
            fdfq.T  = [-flipud(fdfq.T(2:end,:));fdfq.T];
        end
    end
else
    if min(fdfq.Iq,[],'all')>=0
        fdfq.Id = [fdfq.Id(2:end,:);fdfq.Id];
        fdfq.Iq = [-flipud(fdfq.Iq(2:end,:));fdfq.Iq];
        fdfq.Fd = [flipud(fdfq.Fd(2:end,:));fdfq.Fd];
        fdfq.Fq = [-flipud(fdfq.Fq(2:end,:));fdfq.Fq];
        fdfq.T  = [-flipud(fdfq.T(2:end,:));fdfq.T];
    end
end

fd0 = interp2(fdfq.Id,fdfq.Iq,fdfq.Fd,id0,iq0);
fq0 = interp2(fdfq.Id,fdfq.Iq,fdfq.Fq,id0,iq0);
T0  = interp2(fdfq.Id,fdfq.Iq,fdfq.T,id0,iq0);


time = linspace(0,nCycle*t1,nPoints*nCycle+1);

idVect    = nan(size(time));
iqVect    = nan(size(time));
fdVect    = nan(size(time));
fqVect    = nan(size(time));
TVect     = nan(size(time));
TVectExt  = nan(size(time));
intVect   = ones(size(time));

idVect(1)    = id0;
iqVect(1)    = iq0;
fdVect(1)    = fd0;
fqVect(1)    = fq0;
TVect(1)     = T0;
TVectExt(1)  = T0;
intVect(1)   = 1;

% steady-state solution
Idq = fdfq.Id+j*fdfq.Iq;
Fdq = fdfq.Fd+j*fdfq.Fq;
Vdq = Rs*Idq+j*n*pi/30*p*Fdq;
c = contourc(fdfq.Id(1,:),fdfq.Iq(:,1),real(Vdq),[0 0]);
idTmp = c(1,2:end);
iqTmp = c(2,2:end);
iiTmp = 1:1:numel(idTmp);
VqTmp = interp2(fdfq.Id,fdfq.Iq,imag(Vdq),idTmp,iqTmp);
index = interp1(VqTmp,iiTmp,0);

idInf = interp1(iiTmp,idTmp,index);
iqInf = interp1(iiTmp,iqTmp,index);
fdInf = interp2(fdfq.Id,fdfq.Iq,fdfq.Fd,idInf,iqInf);
fqInf = interp2(fdfq.Id,fdfq.Iq,fdfq.Fq,idInf,iqInf);
TInf  = interp2(fdfq.Id,fdfq.Iq,fdfq.T,idInf,iqInf);

%% figures
iStr=num2str(abs(id0+j*iq0),3);
iStr=strrep(iStr,'.','A');
if ~contains(iStr,'A')
    iStr=[iStr 'A'];
end
gammaStr=num2str(atan2(iq0,id0)*180/pi,4);
gammaStr=strrep(gammaStr,'.','d');
if ~contains(gammaStr,'d')
    gammaStr=[gammaStr 'd'];
end
resFolder = [motName '_results\MMM results\' 'TransientSC - ' iStr '_' gammaStr '_' int2str(n) 'rpm' '_' int2str(motorModel.data.tempPM) 'degPM_' int2str(motorModel.data.tempCu) 'degCu' '\'];


hfig(1) = figure();
figSetting();
xlim([min(time) max(time)]);
hax(1) = gca;
xlabel('$t$ [s]')
ylabel('[A]')
hleg(1) = legend(hax(1),'show','Location','southeast');
set(gcf,'FileName',[pathname resFolder 'CurrentVStime.fig'])
plot(hax(1),[time(1) time(end)],idInf*[1 1],'--b','DisplayName','$i_d$ - steady-state')
plot(hax(1),[time(1) time(end)],iqInf*[1 1],'--r','DisplayName','$i_q$ - steady-state')
plot(hax(1),[time(1) time(end)],abs(idInf+j*iqInf)*[1 1],'--k','DisplayName','$|i_{dq}|$ - steady-state')
hp1(1) = plot(hax(1),time,idVect.*intVect,':b','DisplayName','$i_d$','HandleVisibility','off');
hp1(2) = plot(hax(1),time,iqVect.*intVect,':r','DisplayName','$i_q$','HandleVisibility','off');
hp1(3) = plot(hax(1),time,abs(idVect+j*iqVect).*intVect,':k','DisplayName','$|i_{dq}|$','HandleVisibility','off');
hp1(4) = plot(hax(1),time,idVect,'-b','DisplayName','$i_d$');
hp1(5) = plot(hax(1),time,iqVect,'-r','DisplayName','$i_q$');
hp1(6) = plot(hax(1),time,abs(idVect+j*iqVect),'-k','DisplayName','$|i_{dq}|$');


hfig(2) = figure();
figSetting();
set(gca,...
    'DataAspectRatio',[1 1 1],...
    'GridAlpha',1,...
    'Layer','top');
xlabel('$i_d$ [A]')
ylabel('$i_q$ [A]')
title('Flux linkage during short circuit and current trajectory')
set(gcf,'FileName',[pathname resFolder 'idiqTrajectory.fig'])
[c,h] = contourf(fdfq.Id,fdfq.Iq,abs(fdfq.Fd+j*fdfq.Fq),'DisplayName','$|\lambda_{dq}|$');
clabel(c,h);
tmp = get(gca,'CLim');
set(gca,'CLim',tmp,'CLimMode','manual')
[c,h] = contour(fdfq.Id,fdfq.Iq,abs(fdfq.Fd+j*fdfq.Fq),abs(fd0+j*fq0)*[1 1],'-k','LineWidth',1.5,'DisplayName','$|\lambda_{ini}|$');
[c,h] = contour(fdfq.Id,fdfq.Iq,abs(fdfq.Id+j*fdfq.Iq),abs(id0+j*iq0)*[1 1],'--k','LineWidth',1.5,'DisplayName','$|i_{ini}|$');
hax(2) = gca;
hp2(1) = plot(hax(2),idVect,iqVect,'-r','DisplayName','Transient');
plot(hax(2),idInf,iqInf,'r*','DisplayName','Steady-state')
hleg(2) = legend(hax(2),'show','Location','northeast');

hfig(3) = figure();
figSetting();
xlim([min(time) max(time)]);
hax(3) = gca;
xlabel('$t$ [s]')
ylabel('$T$ [Nm]')
hleg(3) = legend(hax(3),'show','Location','southeast');
set(gcf,'FileName',[pathname resFolder 'TorqueVStime.fig'])
plot(hax(3),[time(1) time(end)],TInf*[1 1],'--b','DisplayName','Steady-state')
hp3(1) = plot(hax(3),time,TVect,'-b','DisplayName','Interpolated');
hp3(2) = plot(hax(3),time,TVectExt,':b','DisplayName','Extrapolated');

hfig(4) = figure();
figSetting();
xlim([min(time) max(time)]);
hax(4) = gca;
xlabel('$t$ [s]')
ylabel('$\lambda$ [Vs]')
hleg(4) = legend(hax(4),'show','Location','southeast');
set(gcf,'FileName',[pathname resFolder 'FluxVStime.fig'])
plot(hax(4),[time(1) time(end)],fdInf*[1 1],'--b','DisplayName','$\lambda_d$ - steady-state')
plot(hax(4),[time(1) time(end)],fqInf*[1 1],'--r','DisplayName','$\lambda_q$ - steady-state')
plot(hax(4),[time(1) time(end)],abs(fdInf+j*fqInf)*[1 1],'--k','DisplayName','$|\lambda_{dq}|$ - steady-state')
hp4(1) = plot(hax(4),time,fdVect.*intVect,':b','DisplayName','$\lambda_d$','HandleVisibility','off');
hp4(2) = plot(hax(4),time,fqVect.*intVect,':r','DisplayName','$\lambda_q$','HandleVisibility','off');
hp4(3) = plot(hax(4),time,abs(fdVect+j*fqVect).*intVect,':k','DisplayName','$|\lambda_{dq}|$','HandleVisibility','off');
hp4(4) = plot(hax(4),time,fdVect,'-b','DisplayName','$\lambda_d$');
hp4(5) = plot(hax(4),time,fqVect,'-r','DisplayName','$\lambda_q$');
hp4(6) = plot(hax(4),time,abs(fdVect+j*fqVect),'-k','DisplayName','$|\lambda_{dq}|$');

data.Fd = fdfq.Fd(:);
data.Fq = fdfq.Fq(:);
data.Id = fdfq.Id(:);
data.Iq = fdfq.Iq(:);
data.T  = fdfq.T(:);

filt = ones(1,length(data.Id));
filt(isnan(data.Fd)) = 0;
filt(isnan(data.Fq)) = 0;

data.Fd(filt==0) = [];
data.Fq(filt==0) = [];
data.Id(filt==0) = [];
data.Iq(filt==0) = [];
data.T(filt==0)  = [];

fInt.Id = scatteredInterpolant(data.Fd,data.Fq,data.Id,'linear','none');
fInt.Iq = scatteredInterpolant(data.Fd,data.Fq,data.Iq,'linear','none');
fInt.T  = scatteredInterpolant(data.Fd,data.Fq,data.T,'linear','none');
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
xlabel('$\lambda_d$ [A]')
ylabel('$\lambda_q$ [A]')
title('Current during short circuit and flux linkage trajectory')
set(gcf,'FileName',[pathname resFolder 'fdfqTrajectory.fig'])
[c,h] = contourf(FdPlot,FqPlot,abs(IdPlot+j*IqPlot),'DisplayName','$|i_{dq}|$');
clabel(c,h);
[c,h] = contour(FdPlot,FqPlot,abs(IdPlot+j*IqPlot),abs(id0+j*iq0)*[1 1],'-k','LineWidth',1.5,'DisplayName','$|i_{ini}|$');
tmp = get(gca,'CLim');
set(gca,'CLim',tmp,'CLimMode','manual')
[c,h] = contour(FdPlot,FqPlot,abs(FdPlot+j*FqPlot),abs(fd0+j*fq0)*[1 1],'--k','LineWidth',1.5,'DisplayName','$|\lambda_{ini}|$');
hax(5) = gca;
hp5(1) = plot(hax(5),fdVect,fqVect,'-r','DisplayName','Transient');
plot(hax(5),fdInf,fqInf,'r*','DisplayName','Steady-state')
hleg(5) = legend(hax(5),'show','Location','northeast');

plot(hax(2),idVect(1),iqVect(1),'ro','HandleVisibility','off');
plot(hax(5),fdVect(1),fqVect(1),'ro','HandleVisibility','off');

dt = time(2);
fInt.Id = scatteredInterpolant(data.Fd,data.Fq,data.Id,'linear','linear');
fInt.Iq = scatteredInterpolant(data.Fd,data.Fq,data.Iq,'linear','linear');
fInt.T  = scatteredInterpolant(data.Fd,data.Fq,data.T,'linear','linear');

for ii=1:length(hleg)
    switch ii
        case {3,4,5}
            set(hleg(ii),'NumColumns',2,'Location','southoutside');
        case {1,2}
            set(hleg(ii),'NumColumns',2,'Location','south');
    end
end

for ii=1:length(hfig)
    tmp = get(hfig(ii),'FileName');
    [~,name,~] = fileparts(tmp);
    set(hfig(ii),'Name',name);
end

disp('Transient short-circuit computation...')
fprintf(' %06.2f%%',0)

% tic
%% transient solution
for tt=2:length(time)
    dFd = (+w*fqVect(tt-1)-Rs*idVect(tt-1))*dt;
    dFq = (-w*fdVect(tt-1)-Rs*iqVect(tt-1))*dt;
    fdVect(tt)    = fdVect(tt-1)+dFd;
    fqVect(tt)    = fqVect(tt-1)+dFq;
    idVect(tt)    = fInt.Id(fdVect(tt),fqVect(tt));
    iqVect(tt)    = fInt.Iq(fdVect(tt),fqVect(tt));
    TVect(tt)     = interp2(fdfq.Id,fdfq.Iq,fdfq.T,idVect(tt),iqVect(tt));
    TVectExt(tt)  = fInt.T(fdVect(tt),fqVect(tt));
    intVect = ~isnan(TVect);
    intVect = double(intVect);
    intVect(intVect==0) = NaN;

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
        
        set(hp3(1),'YData',TVect);
        set(hp3(2),'YData',TVectExt);
        
        set(hp4(1),'YData',fdVect);
        set(hp4(2),'YData',fqVect);
        set(hp4(3),'YData',abs(fdVect+j*fqVect));
        set(hp4(4),'YData',fdVect.*intVect);
        set(hp4(5),'YData',fqVect.*intVect);
        set(hp4(6),'YData',abs(fdVect+j*fqVect).*intVect);
        
        set(hp5(1),'XData',fdVect,'YData',fqVect);
        
        drawnow()
    end
end

%toc

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
    
    set(hp3(1),'YData',TVect);
    set(hp3(2),'YData',TVectExt);
    
    set(hp4(1),'YData',fdVect);
    set(hp4(2),'YData',fqVect);
    set(hp4(3),'YData',abs(fdVect+j*fqVect));
    set(hp4(4),'YData',fdVect.*intVect);
    set(hp4(5),'YData',fqVect.*intVect);
    set(hp4(6),'YData',abs(fdVect+j*fqVect).*intVect);
    
    set(hp5(1),'XData',fdVect,'YData',fqVect);
    
    drawnow()
end



% Output data
SCresults.time    = time;
SCresults.id      = idVect;
SCresults.iq      = iqVect;
SCresults.fd      = fdVect;
SCresults.fq      = fqVect;
SCresults.T       = TVectExt;
SCresults.intFilt = intVect;
SCresults.idInf   = idInf;
SCresults.iqInf   = iqInf;
SCresults.fdInf   = fdInf;
SCresults.fqInf   = fqInf;
SCresults.TInf    = TInf;
SCresults.n       = n;

%% save data and figures

answer = 'No';
answer = questdlg('Save results?','Save','Yes','No',answer);
if strcmp(answer,'Yes')
    if ~exist([pathname resFolder],'dir')
        mkdir([pathname resFolder]);
    end
    
    save([pathname resFolder 'TransientShortCircuitResults.mat'],'SCresults')
    for ii=1:length(hfig)
        savePrintFigure(hfig(ii));
    end
end


