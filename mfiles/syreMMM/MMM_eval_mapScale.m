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

function [mapScale] = MMM_eval_mapScale(motorModel,setup)

% Create the maps (kL,kN) of the scaled motor (ECCE 2022)

if nargin==1
    prompt = {'Min axial length','Max axial length','Min number of turns','Max number of turns','DC link voltage [V]','Peak phase current [Apk]'};
    name   = 'Scaling map setup';
    answer = {
        num2str(floor(motorModel.data.l*0.5))
        num2str(ceil(motorModel.data.l*2))
        num2str(floor(motorModel.data.Ns*0.5))
        num2str(ceil(motorModel.data.Ns*2))
        int2str(motorModel.data.Vdc)
        int2str(motorModel.data.Imax)
        };
    answer = inputdlg(prompt,name,1,answer);
    lMin  = eval(answer{1});
    lMax  = eval(answer{2});
    NsMin = eval(answer{3});
    NsMax = eval(answer{4});
    setup.Vdc   = eval(answer{5});
    setup.Imax  = eval(answer{6});
    setup.lVect  = linspace(lMin,lMax,101);
    setup.NsVect = linspace(NsMin,NsMax,101);
end

if isempty(motorModel.controlTrajectories)
    motorModel.controlTrajectories = MMM_eval_AOA(motorModel,'LUT');
end

MTPA = motorModel.controlTrajectories.MTPA;

Imax = setup.Imax;
Vdc  = setup.Vdc;
Rs0  = motorModel.data.Rs;
l0   = motorModel.data.l;
lend = motorModel.data.lend;
Ns0  = motorModel.data.Ns;
p    = motorModel.data.p;
R    = motorModel.data.R;

if strcmp(motorModel.data.motorType,'PM')
    if ~isempty(motorModel.DemagnetizationLimit)
        Idemag0 = interp1(motorModel.DemagnetizationLimit.tempPM,motorModel.DemagnetizationLimit.Idemag,motorModel.data.tempPM);
    else
        Idemag0 = NaN;
    end
else
    Idemag0 = NaN;
end

% evaluation of the HWC current
fdfq = motorModel.FluxMap_dq;
if strcmp(motorModel.data.axisType,'SR')
    iTmp = unique(fdfq.Iq);
    fTmp = interp2(fdfq.Id,fdfq.Iq,fdfq.Fq,zeros(size(iTmp)),iTmp);
else
    iTmp = unique(fdfq.Id);
    fTmp = interp2(fdfq.Id,fdfq.Iq,fdfq.Fd,iTmp,zeros(size(iTmp)));
    iTmp = -iTmp;
    fTmp = -fTmp;
end



[l,Ns] = meshgrid(setup.lVect,setup.NsVect);

T    = nan(size(l));
n    = nan(size(l));
id   = nan(size(l));
iq   = nan(size(l));
fd   = nan(size(l));
fq   = nan(size(l));
iHWC = nan(size(l));

kL = l/l0;
kN = Ns/Ns0;

Rs = kN.^2.*Rs0.*(kL.*l0./(l0+lend)+lend/(l0+lend));

index = 1:1:numel(l);

for ii=1:length(index)
    id_MTPA = MTPA.id/kN(ii);
    iq_MTPA = MTPA.iq/kN(ii);
    T_MTPA  = MTPA.T*kL(ii);
    fd_MTPA = MTPA.fd*kN(ii)*kL(ii);
    fq_MTPA = MTPA.fq*kN(ii)*kL(ii);

    id(ii) = interp1(abs(id_MTPA+j*iq_MTPA),id_MTPA,Imax);
    iq(ii) = interp1(abs(id_MTPA+j*iq_MTPA),iq_MTPA,Imax);
    fd(ii) = interp1(abs(id_MTPA+j*iq_MTPA),fd_MTPA,Imax);
    fq(ii) = interp1(abs(id_MTPA+j*iq_MTPA),fq_MTPA,Imax);
    w_A    = calcLimitPulsation(id(ii),iq(ii),fd(ii),fq(ii),Rs(ii),Vdc/sqrt(3));

    n(ii) = real(w_A)*30/pi/p;
    T(ii) = interp1(abs(id_MTPA+j*iq_MTPA),T_MTPA,Imax);

    iHWC(ii) = interp1(fTmp*kN(ii)*kL(ii),iTmp/kN(ii),abs(fd(ii)+j*fq(ii)),'linear','extrap');
end

loss = 3/2*Rs.*abs(id+j*iq).^2;
kj = loss./(2*pi*R/1000*l/1000);

Idemag = Idemag0./kN;
fM     = interp2(fdfq.Id,fdfq.Iq,abs(fdfq.Fd+j*fdfq.Fq),0,0).*kN.*kL;
nUGO   = n.*abs(fd+j*fq)./fM;

% Save data in the output structure
mapScale.l          = l;
mapScale.Ns         = Ns;
mapScale.kL         = kL;
mapScale.kN         = kN;
mapScale.Rs         = Rs;
mapScale.T          = T;
mapScale.n          = n;
mapScale.id         = id;
mapScale.iq         = iq;
mapScale.fd         = fd;
mapScale.fq         = fq;
mapScale.loss       = loss;
mapScale.kj         = kj;
mapScale.PF         = sin(atan2(mapScale.iq,mapScale.id)-atan2(mapScale.fq,mapScale.fd));
if isfield(motorModel,'geo')
    mapScale.J = abs(mapScale.id+j*mapScale.iq)/sqrt(2)/(motorModel.geo.Aslot*motorModel.geo.win.kcu).*(Ns/(motorModel.data.p*motorModel.geo.q));
else
    mapScale.J = NaN;
end
mapScale.Idemag     = Idemag;
mapScale.iHWC       = iHWC;
mapScale.fM         = fM;
mapScale.nUGO       = nUGO;

mapScale.motorModel = motorModel;

% figure

motName  = motorModel.data.motorName;
pathname = motorModel.data.pathname;

resFolder = [motName '_results\MMM results\MapScalingLN_' datestr(now,30) '\'];

hfig(1) = figure();
figSetting()
xlabel('$L$ [mm]')
ylabel('$N_s$')
colors = get(gca,'ColorOrder');
contour(mapScale.l,mapScale.Ns,mapScale.T,'-','LineColor',colors(2,:),'LineWidth',1.5,'ShowText','on','DisplayName','$T$ [Nm]')
contour(mapScale.l,mapScale.Ns,mapScale.PF,'-','LineColor',colors(1,:),'LineWidth',1.5,'ShowText','on','DisplayName','$cos\varphi$')
if ~isnan(mapScale.J(1,1))
    contour(mapScale.l,mapScale.Ns,mapScale.J,'-','LineColor',colors(3,:),'LineWidth',1.5,'ShowText','on','DisplayName','$J$ [Arms/mm$^2$]')
else
    contour(mapScale.l,mapScale.Ns,mapScale.kj/1000,'-','LineColor',colors(3,:),'LineWidth',1.5,'ShowText','on','DisplayName','$k_j$ [kW/m$^2$]')
end
plot(l0,Ns0,'ko','MarkerFaceColor','k','DisplayName','Baseline')
legend('show','Location','northeast');
title(['($L,N_s)$ map - $V_{dc}=' int2str(Vdc) '$ V / $I_{max}=' int2str(Imax) '$ Apk'])
set(hfig(1),'FileName',[pathname resFolder 'mapScaling.fig'],'UserData',mapScale);

hfig(2) = figure();
figSetting()
xlabel('$L$ [mm]')
ylabel('$N_s$')
colors = get(gca,'ColorOrder');
contour(mapScale.l,mapScale.Ns,mapScale.T,'-','LineColor',colors(2,:),'LineWidth',1.5,'ShowText','on','DisplayName','$T$ [Nm]')
contour(mapScale.l,mapScale.Ns,mapScale.n,'-','LineColor',colors(1,:),'LineWidth',1.5,'ShowText','on','DisplayName','$n$ [rpm]')
contour(mapScale.l,mapScale.Ns,mapScale.kj/1000,'-','LineColor',colors(3,:),'LineWidth',1.5,'ShowText','on','DisplayName','$k_j$ [kW/m$^2$]')
plot(l0,Ns0,'ko','MarkerFaceColor','k','DisplayName','Baseline')
legend('show','Location','northeast');
title(['($L,N_s)$ map - $V_{dc}=' int2str(Vdc) '$ V / $I_{max}=' int2str(Imax) '$ Apk'])
set(hfig(2),'FileName',[pathname resFolder 'mapScaling_T_n_kj.fig'],'UserData',mapScale);

% save data and figures

answer = 'No';
answer = questdlg('Save results?','Save','Yes','No',answer);
if strcmp(answer,'Yes')
    if ~exist([pathname resFolder],'dir')
        mkdir([pathname resFolder]);
    end
    
    save([pathname resFolder 'ScalingMapData.mat'],'mapScale','motorModel')
    for ii=1:length(hfig)
        savePrintFigure(hfig(ii));
    end
end

if nargout()==0
    clear mapScale
end
