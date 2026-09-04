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

% Create the maps (kL,kN) of the scaled motor (ECCE 2022, OJIA 2025)

if nargin==1
    prompt = {'Min axial length','Max axial length','Min number of turns','Max number of turns','DC link voltage (V)','Peak phase current (Apk)','Number of poits per axes'};
    name   = 'Scaling map setup';
    answer = {
        num2str(floor(motorModel.data.l*0.5))
        num2str(ceil(motorModel.data.l*2))
        num2str(floor(motorModel.data.Ns*0.5))
        num2str(ceil(motorModel.data.Ns*2))
        int2str(motorModel.data.Vdc)
        int2str(motorModel.data.Imax)
        int2str(51)
        };
    answer = inputdlg(prompt,name,1,answer);
    lMin  = eval(answer{1});
    lMax  = eval(answer{2});
    NsMin = eval(answer{3});
    NsMax = eval(answer{4});
    nPts  = eval(answer{7});
    setup.Vdc   = eval(answer{5});
    setup.Imax  = eval(answer{6});
    setup.lVect  = linspace(lMin,lMax,nPts);
    setup.NsVect = linspace(NsMin,NsMax,nPts);
end

if isempty(motorModel.controlTrajectories)
    motorModel.controlTrajectories = MMM_eval_AOA(motorModel,'LUT');
end

% MTPA = motorModel.controlTrajectories.MTPA;

Imax = setup.Imax;
Vdc  = setup.Vdc;
Rs0  = motorModel.data.Rs;
l0   = motorModel.data.l;
lend = motorModel.data.lend;
Ns0  = motorModel.data.Ns;
p    = motorModel.data.p;
R    = motorModel.data.R;
n3ph = motorModel.data.n3phase;

motorModel.data.Imax = Imax;
motorModel.data.Vdc  = Vdc;

if ~isempty(motorModel.DemagnetizationLimit)
    tempPMmax = max(motorModel.DemagnetizationLimit.tempPM);
    tempPMmax = 150;
end


if strcmp(motorModel.data.motorType,'PM')
    if ~isempty(motorModel.DemagnetizationLimit)
        Idemag0 = interp1(motorModel.DemagnetizationLimit.tempPM,motorModel.DemagnetizationLimit.Idemag,motorModel.data.tempPM);
        IdemagTempMax0 = interp1(motorModel.DemagnetizationLimit.tempPM,motorModel.DemagnetizationLimit.Idemag,tempPMmax);
        Idemag20deg0 = interp1(motorModel.DemagnetizationLimit.tempPM,motorModel.DemagnetizationLimit.Idemag,20);
    else
        Idemag0 = NaN;
        IdemagTempMax0 = NaN;
        Idemag20deg0 = NaN;
    end
else
    Idemag0 = NaN;
    IdemagTempMax0 = NaN;
    Idemag20deg0 = NaN;
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

T     = nan(size(l));
n     = nan(size(l));
id    = nan(size(l));
iq    = nan(size(l));
fd    = nan(size(l));
fq    = nan(size(l));
iHWC  = nan(size(l));
ich   = nan(size(l));
Pmax  = nan(size(l));
Pnmax = nan(size(l));

kL = l/l0;
kN = Ns/Ns0;

Rs = kN.^2.*Rs0.*(kL.*l0./(l0+lend)+lend/(l0+lend));

index = 1:1:numel(l);

motorModelTmp = motorModel;

disp('(L,Ns) map evaluation in progress...')
fprintf(' %06.2f%%',0)

for ii=1:length(index)

    motorModelTmp.FluxMap_dq.Id = motorModel.FluxMap_dq.Id/kN(ii);
    motorModelTmp.FluxMap_dq.Iq = motorModel.FluxMap_dq.Iq/kN(ii);
    motorModelTmp.FluxMap_dq.Fd = motorModel.FluxMap_dq.Fd*kN(ii)*kL(ii);
    motorModelTmp.FluxMap_dq.Fq = motorModel.FluxMap_dq.Fq*kN(ii)*kL(ii);
    motorModelTmp.FluxMap_dq.T  = motorModel.FluxMap_dq.T*kL(ii);

    motorModelTmp.controlTrajectories.MTPA.id = motorModel.controlTrajectories.MTPA.id/kN(ii);
    motorModelTmp.controlTrajectories.MTPA.iq = motorModel.controlTrajectories.MTPA.iq/kN(ii);
    motorModelTmp.controlTrajectories.MTPA.fd = motorModel.controlTrajectories.MTPA.fd*kN(ii)*kL(ii);
    motorModelTmp.controlTrajectories.MTPA.fq = motorModel.controlTrajectories.MTPA.fq*kN(ii)*kL(ii);
    motorModelTmp.controlTrajectories.MTPA.T  = motorModel.controlTrajectories.MTPA.T*kL(ii);

    motorModelTmp.controlTrajectories.MTPV.id = motorModel.controlTrajectories.MTPV.id/kN(ii);
    motorModelTmp.controlTrajectories.MTPV.iq = motorModel.controlTrajectories.MTPV.iq/kN(ii);
    motorModelTmp.controlTrajectories.MTPV.fd = motorModel.controlTrajectories.MTPV.fd*kN(ii)*kL(ii);
    motorModelTmp.controlTrajectories.MTPV.fq = motorModel.controlTrajectories.MTPV.fq*kN(ii)*kL(ii);
    motorModelTmp.controlTrajectories.MTPV.T  = motorModel.controlTrajectories.MTPV.T*kL(ii);


    id_MTPA = motorModelTmp.controlTrajectories.MTPA.id;
    iq_MTPA = motorModelTmp.controlTrajectories.MTPA.iq;
    T_MTPA  = motorModelTmp.controlTrajectories.MTPA.T;
    fd_MTPA = motorModelTmp.controlTrajectories.MTPA.fd;
    fq_MTPA = motorModelTmp.controlTrajectories.MTPA.fq;

    id(ii) = interp1(abs(id_MTPA+j*iq_MTPA),id_MTPA,Imax);
    iq(ii) = interp1(abs(id_MTPA+j*iq_MTPA),iq_MTPA,Imax);
    fd(ii) = interp1(abs(id_MTPA+j*iq_MTPA),fd_MTPA,Imax);
    fq(ii) = interp1(abs(id_MTPA+j*iq_MTPA),fq_MTPA,Imax);
    % w_A    = calcLimitPulsation(id(ii),iq(ii),fd(ii),fq(ii),Rs(ii),Vdc/sqrt(3));

    % n(ii) = real(w_A)*30/pi/p;
    % T(ii) = interp1(abs(id_MTPA+j*iq_MTPA),T_MTPA,Imax);

    [PlimTmp] = OpLimEval(motorModelTmp,Imax,Vdc);
    Pmax(ii) = max(PlimTmp.P);
    Pnmax(ii) = PlimTmp.P(end);
    n(ii)     = PlimTmp.n_A;
    T(ii)     = PlimTmp.T_A;

    iHWC(ii) = interp1(fTmp*kN(ii)*kL(ii),iTmp/kN(ii),abs(fd(ii)+j*fq(ii)),'linear','extrap');
    ich(ii)  = interp1(fTmp*kN(ii)*kL(ii),iTmp/kN(ii),0,'linear','extrap');

    fprintf('\b\b\b\b\b\b\b')
    fprintf('%06.2f%%',ii/numel(index)*100)
end

disp('(L,Ns) map evaluated!')

loss = 3/2*Rs.*abs(id+j*iq).^2*n3ph;
kj = loss./(2*pi*R/1000*l/1000);

Idemag = Idemag0./kN;
IdemagTempMax = IdemagTempMax0./kN;
Idemag20deg = Idemag20deg0./kN;
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
mapScale.Idemag        = Idemag;
mapScale.IdemagTempMax = IdemagTempMax;
mapScale.Idemag20deg   = Idemag20deg;
mapScale.iHWC          = iHWC;
mapScale.ich           = ich;
mapScale.fM            = fM;
mapScale.nUGO          = nUGO;
mapScale.Pbase         = mapScale.T.*mapScale.n*pi/30;
mapScale.Pmax          = Pmax;
mapScale.Pnmax         = Pnmax;

mapScale.motorModel = motorModel;

% figure

motName  = motorModel.data.motorName;
pathname = motorModel.data.pathname;

resFolder = checkPathSyntax([motName '_results\MMM results\MapScalingLN_' datestr(now,30) '\']);

hfig(1) = figure();
figSetting()
xlabel('$L$ (mm)')
ylabel('$N_s$')
colors = get(gca,'ColorOrder');
contour(mapScale.l,mapScale.Ns,mapScale.T,'-','LineColor',colors(2,:),'LineWidth',1.5,'ShowText','on','DisplayName','$T$ (Nm)')
contour(mapScale.l,mapScale.Ns,mapScale.n,'-','LineColor',colors(1,:),'LineWidth',1.5,'ShowText','on','DisplayName','$n$ (rpm)')
% if ~isnan(mapScale.J(1,1))
%     contour(mapScale.l,mapScale.Ns,mapScale.J,'-','LineColor',colors(3,:),'LineWidth',1.5,'ShowText','on','DisplayName','$J$ (Arms/mm$^2$)')
% else
%     contour(mapScale.l,mapScale.Ns,mapScale.kj/1000,'-','LineColor',colors(3,:),'LineWidth',1.5,'ShowText','on','DisplayName','$k_j$ (kW/m$^2$)')
% end
contour(mapScale.l,mapScale.Ns,mapScale.Pbase/1000,'-','LineColor',colors(3,:),'LineWidth',1.5,'ShowText','on','DisplayName','$P_{base}$ (kW)')
plot(l0,Ns0,'ko','MarkerFaceColor','k','DisplayName','Baseline')
legend('show','Location','northeast');
title(['($L,N_s)$ map - $V_{dc}=' int2str(Vdc) '$ V / $I_{max}=' int2str(Imax) '$ Apk'])
set(hfig(1),'FileName',[pathname resFolder 'mapScaling.fig'],'UserData',mapScale);

hfig(2) = figure();
figSetting()
xlabel('$L$ (mm)')
ylabel('$N_s$')
colors = get(gca,'ColorOrder');
contour(mapScale.l,mapScale.Ns,mapScale.T,'-','LineColor',colors(2,:),'LineWidth',1.5,'ShowText','on','DisplayName','$T_{max}$ (Nm)')
contour(mapScale.l,mapScale.Ns,mapScale.n,'-','LineColor',colors(1,:),'LineWidth',1.5,'ShowText','on','DisplayName','$n_{base}$ (rpm)')
contour(mapScale.l,mapScale.Ns,mapScale.kj/1000,'-','LineColor',colors(3,:),'LineWidth',1.5,'ShowText','on','DisplayName','$k_j$ (kW/m$^2$)')
plot(l0,Ns0,'ko','MarkerFaceColor','k','DisplayName','Baseline')
legend('show','Location','northeast');
title(['($L,N_s)$ map - $V_{dc}=' int2str(Vdc) '$ V / $I_{max}=' int2str(Imax) '$ Apk'])
set(hfig(2),'FileName',[pathname resFolder 'mapScaling_T_n_kj.fig'],'UserData',mapScale);

hfig(3) = figure();
figSetting()
xlabel('$L$ (mm)')
ylabel('$N_s$')
colors = get(gca,'ColorOrder');
contour(mapScale.l,mapScale.Ns,mapScale.T,'-','LineColor',colors(2,:),'LineWidth',1.5,'ShowText','on','DisplayName','$T_{max}$ (Nm)')
contour(mapScale.l,mapScale.Ns,mapScale.n,'-','LineColor',colors(1,:),'LineWidth',1.5,'ShowText','on','DisplayName','$n_{base}$ (rpm)')
contour(mapScale.l,mapScale.Ns,mapScale.Pmax/1000,'-','LineColor',colors(3,:),'LineWidth',1.5,'ShowText','on','DisplayName','$P_{max}$ (kW)')
plot(l0,Ns0,'ko','MarkerFaceColor','k','DisplayName','Baseline')
legend('show','Location','northeast');
title(['($L,N_s)$ map - $V_{dc}=' int2str(Vdc) '$ V / $I_{max}=' int2str(Imax) '$ Apk'])
set(hfig(3),'FileName',[pathname resFolder 'mapScaling_T_n_Pmax.fig'],'UserData',mapScale);

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
