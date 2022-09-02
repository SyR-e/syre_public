% Copyright 2016
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

function [dataSet,flagS] = syrmDesign(dataSet)
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

clc,
[~, ~, geo, per, mat] = data0(dataSet);

% Design equations
switch dataSet.TypeOfRotor
    case 'SPM'
        map = syrmDesign_SPM(dataSet);
    case 'Vtype'
        map = syrmDesign_Vtype(dataSet);
    otherwise
        map = syrmDesign_SyR(dataSet);
end

% FEAfix
if dataSet.FEAfixN==0
    map.kd   = ones(size(map.xx));
    map.kq   = ones(size(map.xx));
    map.km   = ones(size(map.xx));
    map.k0   = ones(size(map.xx));
    map.kg   = ones(size(map.xx));
    map.dg   = zeros(size(map.xx));
    map.kmPM = zeros(size(map.xx));
    map.xRaw = [];
    map.bRaw = [];
else
    [FEAfixOut]=FEAfix(dataSet,geo,map);
    map.kd   = FEAfixOut.kd;
    map.kq   = FEAfixOut.kq;
    map.km   = FEAfixOut.km;
    map.k0   = FEAfixOut.k0;
    map.kg   = FEAfixOut.kg;
    map.dg   = FEAfixOut.dg;
    map.kmPM = FEAfixOut.kmPM;
    map.xRaw = FEAfixOut.xRaw;
    map.bRaw = FEAfixOut.bRaw;
end

if strcmp(dataSet.TypeOfRotor,'SPM')
    map.fd = map.fd.*map.kd;
    map.fq = map.fq.*map.kq;
    %map.gamma = map.gamma.*map.kg;
    map.gamma = map.gamma+map.dg;
    map.id = map.iAmp.*cos(map.gamma*pi/180);
    map.iq = map.iAmp.*sin(map.gamma*pi/180);
    map.T  = 3/2*geo.p*(map.fd.*map.iq-map.fq.*map.id);
    map.PF = abs(sin(atan(map.iq./map.id)-atan(map.fq./map.fd)));
elseif strcmp(dataSet.TypeOfRotor,'Vtype')
    map.fd  = map.fM.*map.km+(map.fd-map.fM).*map.kd;
    map.fq  = map.fq.*map.kq;
    map.fM  = map.fM.*map.km;
    %map.gamma = map.gamma.*map.kg;
    map.gamma = map.gamma+map.dg;
    map.id = map.iAmp.*cos(map.gamma*pi/180);
    map.iq = map.iAmp.*sin(map.gamma*pi/180);
    map.T   = 3/2*geo.p*(map.fd.*map.iq-map.fq.*map.id);
    map.ich = map.ich.*map.km./map.k0;
else
    map.fd = map.fd.*map.kd;
    map.fq = (map.fq+map.fM).*map.kq-map.fM.*map.km;
%     map.fq = map.fq.*map.kq;
    % map.gamma = map.gamma.*map.kg;
    map.gamma = map.gamma+map.dg;
    map.id = map.iAmp.*cos(map.gamma*pi/180);
    map.iq = map.iAmp.*sin(map.gamma*pi/180);
    map.T  = 3/2*geo.p*(map.fd.*map.iq-map.fq.*map.id);
    map.PF = abs(sin(atan(map.iq./map.id)-atan(map.fq./map.fd)));
    map.fM = map.fM.*map.km;
    map.mPM = map.mPM.*map.kmPM;
end

% update flux/N
map.NsI0    = map.i0*dataSet.TurnsInSeries;
map.F0_Ns   = abs(map.fd+j*map.fq)/dataSet.TurnsInSeries;

% Output figure
hfig=figure();
figSetting(15,10)
[c, h] = contour(map.xx,map.bb,map.T,'Color','r','LineWidth',1,'DisplayName','$T$ [Nm]');
clabel(c,h);
if ~strcmp(dataSet.TypeOfRotor,'Vtype')
    [c, h] = contour(map.xx,map.bb,map.PF,0.4:0.02:0.96,'Color','b','LineWidth',1,'DisplayName','$cos \varphi$');
    clabel(c,h);
else
    [c, h] = contour(map.xx,map.bb,map.ich./map.i0,[0:0.1:0.9 1 1.2:0.4:2.8],'Color','b','LineWidth',1,'DisplayName','$\frac{i_{ch}}{i_0}$');
    clabel(c,h);
    [c, h] = contour(map.xx,map.bb,map.ich./map.i0,[1 1],'Color','b','LineWidth',2,'HandleVisibility','off');
end
if ~isempty(map.xRaw)
    plot(map.xRaw,map.bRaw,'Color',[0 0.5 0],'LineStyle','none','Marker','o','MarkerFaceColor',[0 0.5 0],'DisplayName','FEAfix')
end
plot(map.xx(isnan(map.T)),map.bb(isnan(map.T)),'rx','DisplayName','unfeasible','MarkerSize',8)
xlabel('$x$ - rotor / stator split');
if strcmp(dataSet.TypeOfRotor,'SPM')
    ylabel('$l_m/g$ - p.u. magnet size')
elseif strcmp(dataSet.TypeOfRotor,'Vtype')
    ylabel('$hc/g$ - p.u. magnet size')
else
    ylabel('$b$ - p.u. magnetic loading');
end
legend('show','Location','NorthEast')
title('torque and PF tradeoff')

set(hfig,'UserData',map);

% Machine selection
button = questdlg('pick up a machine?','SELECT','Yes','No','Yes');

while isequal(button,'Yes')
    
    figure(hfig)
    [geo.x,geo.b] = ginput(1);
    if strcmp(dataSet.TypeOfRotor,'SPM')||strcmp(dataSet.TypeOfRotor,'Vtype')
        setup=inputdlg({'x','lm/g'},'(x,lm/g) values',1,{num2str(geo.x,3),num2str(geo.b,3)});
    else
        setup=inputdlg({'x','b'},'(x,b) values',1,{num2str(geo.x,3),num2str(geo.b,3)});
    end
    if ~isempty(setup)
        geo.x=eval(setup{1});
        geo.b=eval(setup{2});
    end
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(['x = ' num2str(geo.x) ';']);
    if (strcmp(dataSet.TypeOfRotor,'SPM')||strcmp(dataSet.TypeOfRotor,'Vtype'))
        disp(['lm_g = ' num2str(geo.b) ';']);
    else
        disp(['b = ' num2str(geo.b) ';']);
    end
    disp(['Torque = ' num2str(interp2(map.xx,map.bb,map.T,geo.x,geo.b)) ' Nm;']);
    disp(['PwrFac = ' num2str(interp2(map.xx,map.bb,map.PF,geo.x,geo.b))]);
    disp(['current = ' num2str(interp2(map.xx,map.bb,map.iAmp,geo.x,geo.b)) ' A;']);
    disp(['gamma = ' num2str(interp2(map.xx,map.bb,map.gamma,geo.x,geo.b)) ' deg;']);
    disp(['CurrDens = ' num2str(interp2(map.xx,map.bb,map.J,geo.x,geo.b)) ' A/mm2']);
    disp(['EltLoading = ' num2str(interp2(map.xx,map.bb,map.A,geo.x,geo.b)) ' A/mm']);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    
    % Export to Syre GUI
    geo.r = geo.x*geo.R;                                % rotor radius [mm]
    geo.wt = interp2(map.xx,map.bb,map.wt,geo.x,geo.b); % tooth width
    geo.wt = round(geo.wt*100)/100;
    geo.lt=interp2(map.xx,map.bb,map.lt,geo.x,geo.b);	% slot length
    geo.lt=round(geo.lt*100)/100;
    if strcmp(dataSet.TypeOfRotor,'SPM')
        dataSet.ThicknessOfPM = round(geo.b*geo.g*100)/100;
    elseif strcmp(dataSet.TypeOfRotor,'Vtype')
        hc_pu = interp2(map.xx,map.bb,map.hc_pu,geo.x,geo.b);
        dataSet.HCpu = round(hc_pu*100)/100;
        beta = interp2(map.xx,map.bb,map.beta,geo.x,geo.b);
        geo.betaPMshape = (1-beta*2/pi)*geo.p/(geo.p-1);
        dataSet.betaPMshape = round(geo.betaPMshape*100)/100;
        geo.Ar=interp2(map.xx,map.bb,map.Ar,geo.x,geo.b);                             % shaft radius [mm]
        geo.Ar=round(geo.Ar*100)/100;
        geo.dalpha_pu = interp2(map.xx,map.bb,map.alpha_pu,geo.x,geo.b);
        dataSet.ALPHApu = round(geo.dalpha_pu*100)/100;
        dataSet.PMdim=-1;
    else
        geo.x0=geo.R*geo.x /cos(pi/2/geo.p);
        geo.Ar=interp2(map.xx,map.bb,map.Ar,geo.x,geo.b);                             % shaft radius [mm]
        geo.Ar=round(geo.Ar*100)/100;
        geo.la=interp2(map.xx,map.bb,map.la,geo.x,geo.b);	% total insulation
        
        % hc evaluation - flux barriers design
        geo.alpha=cumsum(geo.dalpha);
        switch dataSet.syrmDesignFlag.hc
            case 0 % hc = cost
                disp('flux barrier design: hc = cost')
            case 1 % pb = cost
                disp('flux barrier design: pbk = sk/hc = cost')
            case 2 % min Lfq
                disp('flux barrier design: hc/(df*sk^0.5) = cost')
        end
        
        switch dataSet.syrmDesignFlag.dx
            case 0 % dx=0
                disp('flux carrier design: dx=0')
            case 1 % constant iron
                disp('flux carrier design: Fe = cost')
            case 2 % iron proportional to first harmonic flux
                disp('flux carrier design: Fe proportional to first harmonic flux')
            case 3 % iron proportional to flux
                disp('flux carrier design: Fe proportional to flux')
        end
        
        geo.hc_pu = zeros(1,geo.nlay);
        geo.dx = zeros(1,geo.nlay);
        [m,n]=size(map.xx);
        hcTmp=zeros(m,n);
        dxTmp=zeros(m,n);
        for ii=1:geo.nlay
            for mm=1:m
                for nn=1:n
                    hcTmp(mm,nn)=map.hc_pu{mm,nn}(ii);
                    dxTmp(mm,nn)=map.dx{mm,nn}(ii);
                end
            end
            geo.hc_pu(ii)=interp2(map.xx,map.bb,hcTmp,geo.x,geo.b);
            geo.dx(ii)=interp2(map.xx,map.bb,dxTmp,geo.x,geo.b);
        end
        
        dataSet.HCpu=round(geo.hc_pu*100)/100;
        dataSet.DepthOfBarrier=round(geo.dx*100)/100;

        dataSet.PMdim = -dataSet.PMdimPU./dataSet.PMdimPU;
        dataSet.PMdim(isnan(dataSet.PMdim)) = 0;
    end
    
    % current phase angle
    temp_id = interp2(map.xx,map.bb,map.id,geo.x,geo.b);  % id [A]
    temp_iq = interp2(map.xx,map.bb,map.iq,geo.x,geo.b);  % iq [A]
    dataSet.GammaPP=round(atan2(temp_iq,temp_id)*180/pi*100)/100;
    dataSet.ThermalLoadKj = interp2(map.xx,map.bb,map.kj,geo.x,geo.b);
    dataSet.CurrentDensity = interp2(map.xx,map.bb,map.J,geo.x,geo.b);
    dataSet.AdmiJouleLosses = dataSet.ThermalLoadKj*(2*pi*dataSet.StatorOuterRadius/1000*dataSet.StackLength/1000);

    switch dataSet.syrmDesignFlag.i0
        case 0
            dataSet.AdmiJouleLosses = NaN;
            dataSet.CurrentDensity  = NaN;
        case 1
            dataSet.AdmiJouleLosses = NaN;
            dataSet.ThermalLoadKj = NaN;
    end
    per.Loss = dataSet.AdmiJouleLosses;
    per.kj   = dataSet.ThermalLoadKj;
    per.J    = dataSet.CurrentDensity;

    Aslots = interp2(map.xx,map.bb,map.Aslots,geo.x,geo.b);
    geo.Aslot = Aslots/(6*geo.p*geo.q*geo.win.n3phase);
    geo.lend = interp2(map.xx,map.bb,map.lend,geo.x,geo.b);

    [per] = calc_i0(geo,per,mat);
%     dataSet.AdmiJouleLosses = per.Loss;
%     dataSet.ThermalLoadKj   = per.kj;
%     dataSet.CurrentDensity  = per.J;
    
    % adjourn dataSet
    dataSet.AirGapRadius=round(geo.r*100)/100;
    dataSet.ShaftRadius = round(geo.Ar*100)/100;
    dataSet.ToothLength=geo.lt;
    dataSet.ToothWidth=geo.wt;
    
    button = questdlg('pick up another machine?','SELECT','Yes','No','Yes');
    
    if isequal(button,'No')
        buttonS = questdlg('save the last machine?','SELECT','Yes','No','Yes');
        figure(hfig)
    end
    
end

if ~exist('buttonS')
    buttonS='No';
end

flagS=0;

if isequal(buttonS,'Yes')
    % save new machine
    flagS=1;
    newnamestring = ['x' num2str(geo.x,2) 'b' num2str(geo.b,2)];
    newnamestring(newnamestring=='.') = '';
    dataSet.currentfilename = strrep(dataSet.currentfilename,'.mat',[newnamestring '.mat']);
end

figure(hfig)

% %
% figure(1);

