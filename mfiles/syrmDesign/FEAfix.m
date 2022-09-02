% Copyright 2019
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

function [FEAfixOut] = FEAfix(dataSet,geo,map,debug)
% 
% [kd,kq] = FEAfix(dataSet,geo,map,setup)
% 

if nargin()==4
    eval_type = debug.eval_type;
    gammaDebug = debug.gammaDebug;
    % useful during debug:
    % - eval_type: standard is 'singt', but in debug you can set also 'flxdn'
    % - gammaDebug: if NaN, the (id,iq) reference from the plane is used
    %               if  0, the simulation is carried in no-load condition
    %               if +1, just the d-axis is excited, with the current computed from the plane
    %               if -1, just the q-axis is excited, with the current computed from the plane
else
    eval_type = 'singt';
    gammaDebug = NaN;
end

% dataSet.kPM=1;
% warning('kPM not read!!!')

gammaFix = dataSet.syrmDesignFlag.gf;

FEAfixN = dataSet.FEAfixN;

clc
disp('FEAfix calibration...')
tic

xx = map.xx;
bb = map.bb;

[m,n]=size(xx);

if isfield(map,'fM')
    fM = map.fM;
else
    fM = zeros(m,n);
end

[~, ~, ~, per, mat] = data0(dataSet);
per.BrPP = dataSet.BrPP;

switch FEAfixN
    case 1
        xRaw=mean(mean(xx));
        bRaw=mean(mean(bb));
    case 4
        xRaw=[xx(1,1),xx(1,end),xx(end,1),xx(end,end)];
        bRaw=[bb(1,1),bb(1,end),bb(end,1),bb(end,end)];
    case 5
        xRaw=[xx(1,1),xx(1,end),xx(end,1),xx(end,end),mean(mean(xx))];
        bRaw=[bb(1,1),bb(1,end),bb(end,1),bb(end,end),mean(mean(bb))];
    case 8
        xMax=max(max(xx));
        xMin=min(min(xx));
        bMax=max(max(bb));
        bMin=min(min(bb));
        xMea=0.5*(xMax+xMin);
        bMea=0.5*(bMax+bMin);
        xRaw=[xMin xMin xMax xMax xMin+(xMea-xMin)/2 xMea xMea+(xMea-xMin)/2 xMea];
        bRaw=[bMin bMax bMax bMin bMea bMea+(bMea-bMin)/2 bMea bMin+(bMea-bMin)/2];
    case 12
        xMax=max(max(xx));
        xMin=min(min(xx));
        bMax=max(max(bb));
        bMin=min(min(bb));
        xTmp=linspace(xMin,xMax,3);
        bTmp=linspace(bMin,bMax,4);
        [xRaw,bRaw] = meshgrid(xTmp,bTmp);
        xRaw = xRaw(:)';
        bRaw = bRaw(:)';
    case 16
        xMax=max(max(xx));
        xMin=min(min(xx));
        bMax=max(max(bb));
        bMin=min(min(bb));
        xTmp=linspace(xMin,xMax,4);
        bTmp=linspace(bMin,bMax,4);
        [xRaw,bRaw] = meshgrid(xTmp,bTmp);
        xRaw = xRaw(:)';
        bRaw = bRaw(:)';
    case 1000
        xRaw=reshape(xx,1,numel(xx));
        bRaw=reshape(bb,1,numel(bb));
%         % FEAfix13
%         xRaw=[xMin xMin xMax xMax xMin+(xMea-xMin)/2 xMea xMea+(xMea-xMin)/2 xMin+(xMea-xMin)/2 xMea xMea+(xMea-xMin)/2 xMin+(xMea-xMin)/2 xMea xMea+(xMea-xMin)/2];
%         bRaw=[bMin bMax bMax bMin bMea+(bMea-bMin)/2 bMea+(bMea-bMin)/2 bMea+(bMea-bMin)/2 bMea bMea bMea bMin+(bMea-bMin)/2 bMin+(bMea-bMin)/2 bMin+(bMea-bMin)/2];
%         % FEAfix9
%         xTmp = linspace(xMin,xMax,5);
%         bTmp = linspace(bMin,bMax,5);
%         xTmp = xTmp(2:4);
%         bTmp = bTmp(2:4);
%         [xRaw,bRaw] = meshgrid(xTmp,bTmp);
%         xRaw = xRaw(:)';
%         bRaw = bRaw(:)';
    otherwise
        error('Put a correct number!!!')
end

errorFlag=zeros(size(xRaw));

if strcmp(geo.RotType,'SPM')
    RQnames{1}='hc';
    RQnames{2}='r';
    RQnames{3}='wt';
    RQnames{4}='lt';
    RQnames{5}='gamma';
elseif strcmp(geo.RotType,'Vtype')
    RQnames{1}='dalpha';
    RQnames{2}='hc_pu';
    RQnames{3}='r';
    RQnames{4}='wt';
    RQnames{5}='lt';
    RQnames{6}='betaPMshape';
    RQnames{7}='gamma';
    geo.PMdim = -1;
else
    index=0;
    for ii=1:geo.nlay
        RQnames{index+ii}=['hc_pu(' int2str(ii) ')'];
    end
    index=length(RQnames);
    for ii=1:geo.nlay
        RQnames{index+ii}=['dx(' int2str(ii) ')'];
    end
    index=length(RQnames);
    RQnames{index+1}='r';
    RQnames{index+2}='wt';
    RQnames{index+3}='lt';
    index = index+3;
    for ii=1:numel(geo.PMdim)
        RQnames{index+ii} = ['PMdim(' int2str(ii) ')'];
    end
    index=length(RQnames);
    RQnames{index+1}='gamma';
end

geo.RQnames=RQnames;

RQ     = zeros(length(RQnames),length(xRaw));
fdFEA  = zeros(1,length(xRaw));
fqFEA  = zeros(1,length(xRaw));
fMFEA  = zeros(1,length(xRaw));
f0FEA  = zeros(1,length(xRaw));
gFEA   = zeros(1,length(xRaw));
mPMFEA = zeros(1,length(xRaw));
fdMod  = zeros(1,length(xRaw));
fqMod  = zeros(1,length(xRaw));
fMMod  = zeros(1,length(xRaw));
f0Mod  = zeros(1,length(xRaw));
gMod   = zeros(1,length(xRaw));
mPMMod = zeros(1,length(xRaw));
% i0    = zeros(1,length(xRaw));

if strcmp(eval_type,'flxdn')
    BtFEA = zeros(1,length(xRaw));
    ByFEA = zeros(1,length(xRaw));
    BgFEA = zeros(1,length(xRaw));
end


OBJnames{1} = 'Torque';
geo.OBJnames=OBJnames;
per.objs = [per.min_exp_torque 1 0];

per.nsim_singt      = 6;
per.delta_sim_singt = 60;
per.overload        = 1;
geo.mesh_K          = 5;
geo.mesh_K_MOOA     = 10;
per.BrPP            = dataSet.Br;
per.tempPP          = dataSet.PMtemp;

if dataSet.syrmDesignFlag.i0==0
    per.J    = NaN;
    per.kj   = dataSet.ThermalLoadKj;
    per.Loss = NaN;
else
    per.J    = dataSet.CurrentDensity;
    per.kj   = NaN;
    per.Loss = NaN;
end


for mot=1:length(xRaw)
    %disp([' - Machine design ' int2str(mot) ' of ' int2str(length(xRaw))])
    geo.x = xRaw(mot);
    geo.b = bRaw(mot);

    r = geo.x*geo.R;                                        % rotor radius [mm]
    wt = interp2(xx,bb,map.wt,geo.x,geo.b);                       % tooth width
    wt = round(wt*100)/100;
    lt=interp2(xx,bb,map.lt,geo.x,geo.b);                         % slot length
    lt=round(lt*100)/100;
    %geo.x0=geo.R * geo.x /cos(pi/2/geo.p);
    Ar(mot)=interp2(xx,bb,map.Ar,geo.x,geo.b);                          % shaft radius [mm]
    Ar(mot)=round(Ar(mot)*100)/100;
    %geo.la=interp2(xx,bb,map.la,geo.x,geo.b);                         % total insulation
    
    % current phase angle
    temp_id = interp2(xx,bb,map.id,geo.x,geo.b);         % id [A]
    temp_iq = interp2(xx,bb,map.iq,geo.x,geo.b);         % iq [A]
    %gamma   = round(atan2(temp_iq,temp_id)*180/pi*100)/100;
    gamma = atan2(temp_iq,temp_id)*180/pi;
    
    if ~isnan(gammaDebug)
        switch gammaDebug
            case 0 % no-load
                per.overload = 0;
                gamma = 0;
            case +1
                per.overload = cosd(gamma);
                gamma = 0;
            case +1
                per.overload = sind(gamma);
                gamma = 90;
        end
    end
    
    if (isnan(lt)||(temp_iq==0))
        errorFlag(mot)=1;
    else
        switch geo.RotType
            case 'SPM'
                geo.hc_pu=interp2(xx,bb,map.bb*geo.g,geo.x,geo.b);
                per.i0 = interp2(map.xx,map.bb,map.i0,geo.x,geo.b);
                % fill RQ
                RQ(:,mot)=[geo.hc_pu r wt lt gamma]';
            case 'Vtype'
                geo.dalpha_pu = interp2(map.xx,map.bb,map.alpha_pu,geo.x,geo.b);
                map.beta_pu = (1-map.beta*2/pi)*geo.p/(geo.p-1);
                geo.betaPMshape  = interp2(xx,bb,map.beta_pu,geo.x,geo.b);
                geo.hc_pu        = interp2(xx,bb,map.hc_pu,geo.x,geo.b);
                geo.PMdim = -1;
                per.i0 = interp2(map.xx,map.bb,map.i0,geo.x,geo.b);
                % fill RQ
                RQ(:,mot)=[geo.dalpha_pu geo.hc_pu r wt lt geo.betaPMshape gamma]';
            otherwise
                geo.x0=geo.R * geo.x /cos(pi/2/geo.p);
                geo.la=interp2(xx,bb,map.la,geo.x,geo.b);                         % total insulation
                hcTmp=zeros(m,n);
                dxTmp=zeros(m,n);
                for ii=1:geo.nlay
                    for mm=1:m
                        for nn=1:n
                            hcTmp(mm,nn)=map.hc_pu{mm,nn}(ii);
                            dxTmp(mm,nn)=map.dx{mm,nn}(ii);
                        end
                    end
                    geo.hc_pu(ii)=interp2(xx,bb,hcTmp,geo.x,geo.b);
                    geo.dx(ii)=interp2(xx,bb,dxTmp,geo.x,geo.b);
                end
                % fill RQ
                per.i0 = interp2(map.xx,map.bb,map.i0,geo.x,geo.b);
                PMdimPU = dataSet.PMdimPU./dataSet.PMdimPU;
                PMdimPU(isnan(PMdimPU)) = 0;
                PMdimPU = dataSet.kPM*PMdimPU;
                PMdim = PMdimPU;

                RQ(:,mot)=[geo.hc_pu geo.dx r wt lt PMdim(:)' gamma]';
        end
        
        % model flux linkages
        fdMod(mot) = interp2(xx,bb,map.fd,geo.x,geo.b); % fd [Vs]
        fqMod(mot) = interp2(xx,bb,map.fq,geo.x,geo.b); % fq [Vs]
        gMod(mot)  = gamma; % gamma [degrees]
        if strcmp(geo.RotType,'Vtype')
            f0Mod(mot) = interp2(xx,bb,map.fM-map.Lbase.*(map.Ldpu+map.Lspu).*map.iAmp,geo.x,geo.b);
        else
            f0Mod(mot) = 0;
        end
        if isfield(map,'fM')
            fMMod(mot) = interp2(xx,bb,map.fM,geo.x,geo.b); % fM [Vs]
        end
        if isfield(map,'mPM')
            mPMMod(mot) = interp2(map.xx,map.bb,map.mPM,geo.x,geo.b);
        else
            mPMMod(mot) = 1;
        end
    end
end

geo.Ar = min([Ar,geo.Ar]);

RQ=RQ(:,~errorFlag); % filter the unfeasible machines
xRaw=xRaw(~errorFlag);
bRaw=bRaw(~errorFlag);

save('dataSet','dataSet');

filemot = [dataSet.currentpathname dataSet.currentfilename];
filemot = strrep(filemot,'.mat','.fem');



if ~isempty(RQ)
    ppState=parallelComputingCheck();
    if ppState==0 && FEAfixN==1000
        disp('Parallel computing not enabled. Enabling...')
        parpool();
        ppState=parallelComputingCheck();
        disp(['Parallel computing enabled on ' int2str(ppState) ' workers']);
    end
    if ppState<1
        for mot=1:length(xRaw)
            % disp([' - FEA simulation ' int2str(mot) ' of ' int2str(length(xRaw))])
            FEA=FEAfixSimulation(RQ(:,mot),geo,per,mat,eval_type,filemot,gammaFix);
            fdFEA(mot)  = FEA.fd;
            fqFEA(mot)  = FEA.fq;
            fMFEA(mot)  = FEA.fM;
            f0FEA(mot)  = FEA.f0;
            BgFEA(mot)  = FEA.Bg;
            BtFEA(mot)  = FEA.Bt;
            ByFEA(mot)  = FEA.By;
            gFEA(mot)   = atan2(FEA.iq,FEA.id)*180/pi;
            mPMFEA(mot) = FEA.mPM;
            %mCuFEA(mot) = FEA.mCu;
            disp([' - motor ' int2str(mot) ' of ' int2str(length(xRaw)) ' evaluated with ' int2str(FEA.nFEA) ' FEA'])
        end
    else
        parfor mot=1:length(xRaw)
            % disp([' - FEA simulation ' int2str(mot) ' of ' int2str(length(xRaw))])
            FEA=FEAfixSimulation(RQ(:,mot),geo,per,mat,eval_type,filemot,gammaFix);
            fdFEA(mot)  = FEA.fd;
            fqFEA(mot)  = FEA.fq;
            fMFEA(mot)  = FEA.fM;
            f0FEA(mot)  = FEA.f0;
            BgFEA(mot)  = FEA.Bg;
            BtFEA(mot)  = FEA.Bt;
            ByFEA(mot)  = FEA.By;
            gFEA(mot)   = atan2(FEA.iq,FEA.id)*180/pi;
            mPMFEA(mot) = FEA.mPM;
            %mCuFEA(mot) = FEA.mCu;
            disp([' - motor ' int2str(mot) ' of ' int2str(length(xRaw)) ' evaluated with ' int2str(FEA.nFEA) ' FEA'])
        end
    end
else
    disp('FEAfix motor on the x-b plane are not drawable')
    disp('Please correct the x and b ranges to improve the FEA calibration')
    warning('Calibration error')
end

if strcmp(geo.RotType,'SPM')
    kdRaw = fdFEA./fdMod;
    kqRaw = fqFEA./fqMod;
    kmRaw = zeros(size(fdFEA));
    k0Raw = zeros(size(fdFEA));
elseif strcmp(geo.RotType,'Vtype')
    kmRaw = fMFEA./fMMod;
    k0Raw = (f0FEA-fMFEA)./(f0Mod-fMMod);
    kdRaw = (fdFEA-fMFEA)./(fdMod-fMMod);
    kqRaw = fqFEA./fqMod;
else
    kdRaw   = fdFEA./fdMod;
%     kqRaw = fqFEA./fqMod;
    kqRaw   = (fqFEA+fMFEA)./(fqMod+fMMod);
    kmRaw   = fMFEA./fMMod;
    k0Raw   = zeros(size(fdFEA));
    kmPMRaw = mPMFEA./mPMMod;
end

if gammaFix
    kgRaw = gFEA./gMod;
    dgRaw = gFEA-gMod;
else
    kgRaw = ones(size(kdRaw));
    dgRaw = zeros(size(kdRaw));
end

kdRaw = kdRaw(~errorFlag);
kqRaw = kqRaw(~errorFlag);
kmRaw = kmRaw(~errorFlag);
k0Raw = k0Raw(~errorFlag);
kgRaw = kgRaw(~errorFlag);
dgRaw = dgRaw(~errorFlag);
kmPMRaw = kmPMRaw(~errorFlag);

if FEAfixN==1
    kd = kdRaw*ones(size(xx));
    kq = kqRaw*ones(size(xx));
    km = kmRaw*ones(size(xx));
    k0 = k0Raw*ones(size(xx));
    kg = kgRaw*ones(size(xx));
    dg = dgRaw*ones(size(xx));
    kmPM = kmPMRaw*ones(size(xx));
else
    kd   = scatteredInterpolant(xRaw',bRaw',kdRaw','linear');
    kq   = scatteredInterpolant(xRaw',bRaw',kqRaw','linear');
    km   = scatteredInterpolant(xRaw',bRaw',kmRaw','linear');
    k0   = scatteredInterpolant(xRaw',bRaw',k0Raw','linear');
    kg   = scatteredInterpolant(xRaw',bRaw',kgRaw','linear');
    dg   = scatteredInterpolant(xRaw',bRaw',dgRaw','linear');
    kmPM = scatteredInterpolant(xRaw',bRaw',kmPMRaw','linear');
    kd   = kd(xx,bb);
    kq   = kq(xx,bb);
    km   = km(xx,bb);
    k0   = k0(xx,bb);
    kg   = kg(xx,bb);
    dg   = dg(xx,bb);
    kmPM = kmPM(xx,bb);
end

disp('End of FEAfix calibration')

delete('dataSet.mat')

FEAfixOut.kd   = kd;
FEAfixOut.kq   = kq;
FEAfixOut.km   = km;
FEAfixOut.k0   = k0;
FEAfixOut.kg   = kg;
FEAfixOut.dg   = dg;
FEAfixOut.kmPM = kmPM;
FEAfixOut.xRaw = xRaw;
FEAfixOut.bRaw = bRaw;

if strcmp(eval_type,'flxdn')
    FEAfixOut.Bg = BgFEA;
    FEAfixOut.Bt = BtFEA;
    FEAfixOut.By = ByFEA;
end

if (0)
    figure()
    figSetting(12,12)
    xlabel('$x$')
    ylabel('$b$')
    zlabel('$k_d$')
    view(45,30)
    set(gca,'XLim',[min(min(xx)) max(max(xx))],'YLim',[min(min(bb)) max(max(bb))]);
    surf(xx,bb,kd);
    plot3(xRaw,bRaw,kdRaw,'ko')
    
    figure()
    figSetting(12,12)
    xlabel('$x$')
    ylabel('$b$')
    zlabel('$k_d$')
    view(45,30)
    set(gca,'XLim',[min(min(xx)) max(max(xx))],'YLim',[min(min(bb)) max(max(bb))]);
    surf(xx,bb,kq);
    plot3(xRaw,bRaw,kqRaw,'ko')
    
    keyboard
end

timeEnd=toc();
%clc
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['FEAfix procedure'])
disp(['- number of FEA machines : ' int2str(length(xRaw))])
disp(['- elapsed time           : ' num2str(round(timeEnd,1)) ' s'])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ')

%keyboard

