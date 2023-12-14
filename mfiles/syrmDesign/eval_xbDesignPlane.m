% Copyright 2023
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

function [dataSet,flagS,hfig,map] = eval_xbDesignPlane(dataSet,debug)
%   syrmDesign
%   Script for preliminary design of a Synchonous Reluctance Machine (SyRM)

%   The equations follow the the literature of closed-form design of SyRMs.
%   Main reference are this tutorial course notes:
%   Lipo, T. A., et al. "Synchronous reluctance drives tutorial." IEEE-IAS Annual Meeting. 1994
%   Chapter 3, presented by Prof. A. Vagati, is the one to go and look for
%
%   syrmDesign produces a parametric study, function of x and b
%   one (x,b) or (x, lm/g) combination can be selected from the figure
%   one machine will be saved and visualized in syre

% debug mode controlled with the variable debug:
%  - debug = 0 --> work from GUI: motor selection is implemented (default)
%  - debug = 1 --> lauch outside GUI, motor selection is not implemented

if nargin()==1
    debug = 0;
end

clc,
[~, ~, geo, per, mat] = data0(dataSet);

% Default syrmDesign flag
dataSet.syrmDesignFlag.ks = 1;  % always account for saturation factor
if strcmp(dataSet.FluxBarrierMaterial,'Air')
    dataSet.syrmDesignFlag.ichf     = 0;  % do not compute ich if SyRM
    dataSet.syrmDesignFlag.scf      = 0;  % do not compute iHWC if SyRM
    dataSet.syrmDesignFlag.demag0   = 0;  % do not compute demagnetization if SyRM
    dataSet.syrmDesignFlag.demagHWC = 0;  % do not compute demagnetization if SyRM
end

dataSet.syrmDesignFlag.mech  = 0;   %if 1, structural FEA are run in the FEAfix points and the avg stress in the ribs is retrieved
dataSet.syrmDesignFlag.therm = 0;   %if 1, thermal transient simulations with Motor-CAD are run in the FEAfix points and the max copper temperature is retrived

% Design equations
% switch dataSet.TypeOfRotor
%     case 'SPM'
%         map = syrmDesign_SPM(dataSet);
%     case 'Vtype'
%         map = syrmDesign_Vtype(dataSet);
%     otherwise
%         map = syrmDesign_SyR(dataSet);
% end
map = xbPlane_analyticalDesign(dataSet);

% FEAfix
if dataSet.FEAfixN==0
    map.kd        = ones(size(map.xx));
    map.kq        = ones(size(map.xx));
    map.km        = ones(size(map.xx));
    map.k0        = ones(size(map.xx));
    map.kg        = ones(size(map.xx));
    map.dg        = zeros(size(map.xx));
    map.kmPM      = zeros(size(map.xx));
    map.kich      = zeros(size(map.xx));
    map.kiHWC     = zeros(size(map.xx));
    map.kBmin0    = zeros(size(map.xx));
    map.kdPM0     = zeros(size(map.xx));
    map.kBminHWC  = zeros(size(map.xx));
    map.kdPMHWC   = zeros(size(map.xx));

%     map.kmechrad      = cell(size(map.xx));
%     map.kmechtan      = cell(size(map.xx));
    map.ktempCuMax    = zeros(size(map.xx));
    map.ktempCuMaxAct = zeros(size(map.xx));

    map.xRaw      = [];
    map.bRaw      = [];
else
    [FEAfixOut] = FEAfix(dataSet,geo,map);
    map.kd       = FEAfixOut.kd;
    map.kq       = FEAfixOut.kq;
    map.km       = FEAfixOut.km;
    map.k0       = FEAfixOut.k0;
    map.kg       = FEAfixOut.kg;
    map.dg       = FEAfixOut.dg;
    map.kmPM     = FEAfixOut.kmPM;
    map.kich     = FEAfixOut.kich;
    map.kiHWC    = FEAfixOut.kiHWC;
    map.kBmin0   = FEAfixOut.kBmin0;
    map.kdPM0    = FEAfixOut.kdPM0;
    map.kBminHWC = FEAfixOut.kBminHWC;
    map.kdPMHWC  = FEAfixOut.kdPMHWC;

    map.kmechrad      = FEAfixOut.kmechrad;
    map.kmechtan      = FEAfixOut.kmechtan;
    map.ktempCuMax    = FEAfixOut.kthetaCu;
    map.ktempCuMaxAct = FEAfixOut.kthetaCuAct;

    map.xRaw     = FEAfixOut.xRaw;
    map.bRaw     = FEAfixOut.bRaw;
end

if strcmp(dataSet.TypeOfRotor,'SPM')
    map.fd      = map.fd.*map.kd;
    map.fq      = map.fq.*map.kq;
    map.fM      = map.fM.*map.km;
    map.gamma   = map.gamma+map.dg;
    map.id      = map.iAmp.*cos(map.gamma*pi/180);
    map.iq      = map.iAmp.*sin(map.gamma*pi/180);
    map.T       = 3/2*geo.p*(map.fd.*map.iq-map.fq.*map.id)*geo.win.n3phase;
    map.PF      = abs(sin(atan(map.iq./map.id)-atan(map.fq./map.fd)));
    map.mPM     = map.mPM.*map.kmPM;
    map.ich     = map.ich.*map.kich;
    map.iHWC    = map.iHWC.*map.kiHWC;
    map.Bmin0   = map.Bmin0.*map.kBmin0;
    map.dPM0    = map.dPM0.*map.kdPM0;
    map.BminHWC = map.BminHWC.*map.kBminHWC;
    map.dPMHWC  = map.dPMHWC.*map.kdPMHWC;
elseif strcmp(dataSet.TypeOfRotor,'Vtype')
    map.fd  = map.fM.*map.km+(map.fd-map.fM).*map.kd;
    map.fq  = map.fq.*map.kq;
    map.fM  = map.fM.*map.km;
    map.gamma = map.gamma+map.dg;
    map.id = map.iAmp.*cos(map.gamma*pi/180);
    map.iq = map.iAmp.*sin(map.gamma*pi/180);
    map.T   = 3/2*geo.p*(map.fd.*map.iq-map.fq.*map.id)*geo.win.n3phase;
    map.ich = map.ich.*map.km./map.k0;
else
    map.fd      = map.fd.*map.kd;
    map.fq      = (map.fq+map.fM).*map.kq-map.fM.*map.km;
    map.gamma   = map.gamma+map.dg;
    map.id      = map.iAmp.*cos(map.gamma*pi/180);
    map.iq      = map.iAmp.*sin(map.gamma*pi/180);
    map.T       = 3/2*geo.p*(map.fd.*map.iq-map.fq.*map.id)*geo.win.n3phase;
    map.PF      = abs(sin(atan(map.iq./map.id)-atan(map.fq./map.fd)));
    map.fM      = map.fM.*map.km;
    map.mPM     = map.mPM.*map.kmPM;
    map.ich     = map.ich.*map.kich;
    map.iHWC    = map.iHWC.*map.kiHWC;
    map.Bmin0   = map.Bmin0.*map.kBmin0;
    map.dPM0    = map.dPM0.*map.kdPM0;
    map.BminHWC = map.BminHWC.*map.kBminHWC;
    map.dPMHWC  = map.dPMHWC.*map.kdPMHWC;

    %     map.mechStressRad = map.mechStressRad.*map.kmechrad;
    %     map.mechStressTan = map.mechStressTan.*map.kmechtan;
    [Mraw,Mcol] = size(map.mechStressRad);
    for ii=1:Mraw
        for jj=1:Mcol
            map.mechStressRad{ii,jj}    = map.mechStressRad{ii,jj}.*map.kmechrad{ii,jj};
            map.mechStressTan{ii,jj}    = map.mechStressTan{ii,jj}.*map.kmechtan{ii,jj};
        end
    end
    map.tempCuMax     = (map.dTempCu+per.temphous).*map.ktempCuMax;
    map.tempCuMaxAct  = (map.dTempCu+per.temphous).*map.ktempCuMaxAct;
end

% update flux/N, Nsi0, current loading and UGO factor
map.NsI0    = map.i0*dataSet.TurnsInSeries;
map.F0_Ns   = abs(map.fd+j*map.fq)/dataSet.TurnsInSeries;
map.NsIch   = map.ich*dataSet.TurnsInSeries;
map.NsIHWC  = map.iHWC*dataSet.TurnsInSeries;
map.Ach     = (6*dataSet.TurnsInSeries*map.ich)./(2*pi*dataSet.StatorOuterRadius.*map.xx);
map.AHWC    = (6*dataSet.TurnsInSeries*map.iHWC)./(2*pi*dataSet.StatorOuterRadius.*map.xx);
map.fUGOpu  = abs(map.fd+j*map.fq)./map.fM;

if ~dataSet.syrmDesignFlag.demag0
    map = rmfield(map,'Bmin0');
    map = rmfield(map,'dPM0');
end

if ~dataSet.syrmDesignFlag.demagHWC
    map = rmfield(map,'BminHWC');
    map = rmfield(map,'dPMHWC');
end

if ~dataSet.syrmDesignFlag.mech
    map = rmfield(map,'kmechrad');
    map = rmfield(map,'kmechtan');
    map = rmfield(map,'mechStressRad');
    map = rmfield(map,'mechStressTan');
end

if ~dataSet.syrmDesignFlag.therm
    map = rmfield(map,'ktempCuMax');
    map = rmfield(map,'tempCuMax');
    map = rmfield(map,'ktempCuMaxAct');
    map = rmfield(map,'tempCuMaxAct');
end

% Output figure
hfig=figure();
figSetting(15,10)
contour(map.xx,map.bb,map.T,'Color','r','LineWidth',1,'DisplayName','$T$ [Nm]','ShowText','on');
if ~strcmp(dataSet.TypeOfRotor,'Vtype')
    contour(map.xx,map.bb,map.PF,0.4:0.02:0.96,'Color','b','LineWidth',1,'DisplayName','$cos \varphi$','ShowText','on');
else
    contour(map.xx,map.bb,map.ich./map.i0,[0:0.1:0.9 1 1.2:0.4:2.8],'Color','b','LineWidth',1,'DisplayName','$\frac{i_{ch}}{i_0}$','ShowText','on');
    contour(map.xx,map.bb,map.ich./map.i0,[1 1],'Color','b','LineWidth',2,'HandleVisibility','off');
end
if ~isempty(map.xRaw)
    plot(map.xRaw,map.bRaw,'Color',[0 0.5 0],'LineStyle','none','Marker','o','MarkerFaceColor',[0 0.5 0],'DisplayName','FEAfix')
end
plot(map.xx(isnan(map.T)),map.bb(isnan(map.T)),'rx','DisplayName','unfeasible','MarkerSize',8)
xlabel('$x$ - rotor / stator split');
% if strcmp(dataSet.TypeOfRotor,'SPM')
%     ylabel('$l_m/g$ - p.u. magnet size')
% elseif strcmp(dataSet.TypeOfRotor,'Vtype')
%     ylabel('$hc/g$ - p.u. magnet size')
% else
    ylabel('$b$ - p.u. magnetic loading');
% end
legend('show','Location','NorthEast')
title('torque and PF tradeoff')

set(hfig,'UserData',map);

if ~debug
    button = questdlg('Open in (x,b) Design Plane Explorer?','SELECT','Yes','No','Yes');
    if strcmp(button,'Yes')
        hApp = findall(0,'Name','SyR-e (x,b) Design Plane Explorer');
        if isempty(hApp)
            xbDesignPlaneExplorer(map);
        else
            hApp.RunningAppInstance.map = map;
            hApp.RunningAppInstance.SDE_update;
            figure(hApp)
%             hApp.WindowStyle = 'alwaysontop';
%             hApp.WindowStyle = 'normal';
            clear hApp;
        end
    end
end


flagS = 0;

figure(hfig)
