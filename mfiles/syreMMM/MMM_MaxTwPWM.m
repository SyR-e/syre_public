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

function MMM_MaxTwPWM(motorModel,hax)
%% Input

prompt = {'Speed steps','Torque steps','Speed range','Torque range'};
name   = 'PWM Efficiency Map Setup';
answer = {
    num2str(6)
    num2str(5)
    num2str([300*motorModel.SyreDrive.Converter.fPWM/500/motorModel.data.p motorModel.data.nmax])     %%% speed that corresponds to 500 FEA steps if the fPWM has 10 samples per period
    num2str([30 round(motorModel.data.T0)])
    };
answer = inputdlg(prompt,name,1,answer);
nstep  = eval(answer{1});
Tstep  = eval(answer{2});
nrange = str2num(answer{3});
Trange = str2num(answer{4});


%EXTRACT DATA
%TwMap
fPWM      = motorModel.SyreDrive.Converter.fPWM;
nmin      = nrange(1);
nmax      = nrange(2);
Tmin      = Trange(1);
Tmax      = Trange(2);
% nminSIN   = motorModel.TnSetup.nmin;
% nmaxSIN   = motorModel.TnSetup.nmax;
% TmaxSIN   = motorModel.TnSetup.Tmax;
% TminSIN   = motorModel.TnSetup.Tmin;

% General
p        = motorModel.geo.p;
Imax     = motorModel.data.Imax;
Vdc      = motorModel.data.Vdc;
pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
filename = [motName '.fem'];

%Copper Loss
kac   = motorModel.acLossFactor.k;
fkac  = motorModel.acLossFactor.f;
Tkac  = motorModel.acLossFactor.T;
Tcu   = motorModel.TnSetup.temperature;
Rs0   = motorModel.data.Rs;
temp0 = motorModel.data.tempCu;
l     = motorModel.data.l;
lend  = motorModel.data.lend;


%% Calculate Operating Limits
motorModel.data.nmax = motorModel.data.nmax*1.3;
Plim = OpLimEval(motorModel,Imax,Vdc);
motorModel.data.nmax = motorModel.data.nmax/1.3;

%% Speed chosen according to PWM frequency
if nmin == 0
    nmin = 100;
end
n_vect = nmin:(nmax-nmin)/(nstep-1):nmax;

count1 = floor(60/nmax/p*fPWM);
count2 = floor(60/nmin/p*fPWM);
f_eval = fPWM./linspace(count1,count2,count2-count1);
n_PWM = 60*f_eval/p;

n_eval = zeros(1,nstep);
for ii=1:nstep
    [~,pos] = min(abs(n_PWM-n_vect(ii)));
    n_eval(ii) = n_PWM(pos);
end

%% Torque and Speed grids for PWM effy maps
Tlim1 = Plim.T;
nlim  = Plim.n;

Tlim  = interp1(nlim,Tlim1*0.85,n_eval);    %% T evaluated until 90% of the max value at each speed
% T_grid = (Tlim'*linspace(Tmin/max(Tlim),1,Tstep))';
for ii=1:nstep
    T_grid(:,ii) = (Tlim(ii)*linspace(Tmin/Tlim(ii),1,Tstep))';
end
n_grid = repmat(n_eval,Tstep,1);

%% GUI Axis Plot
if nargin==1
    % no link to an existing axis...create a new axis
    figure()
    figSetting();
    hax = axes('OuterPosition',[0 0 1 1]);
end
cla(hax)
set(hax,'XLim',[0 max(n_eval)]);
set(hax,'YLim',[0 Tmax]);

hgrid  = plot(hax,0,0,'k.','LineWidth',1.5);
hdrive = plot(hax,0,0,'r.','LineWidth',1.5);
hfemm  = plot(hax,0,0,'ko','LineWidth',1.5);

for ii=1:numel(n_grid)
    xdata = [get(hgrid,'XData') n_grid(ii)];
    ydata = [get(hgrid,'YData') T_grid(ii)];
    set(hgrid,'XData',xdata,'YData',ydata);
    drawnow();
end

%% Syredrive
Iwave = cell(Tstep,nstep);

oldpath = pwd;
cd([motorModel.data.pathname motorModel.data.motorName '_ctrl_INST'])

syreDriveSingt.tSim = 0.45;

for ii=1:numel(T_grid)
    %plot
    xdata = [get(hdrive,'XData') n_grid(ii)];
    ydata = [get(hdrive,'YData') T_grid(ii)];
    set(hdrive,'XData',xdata,'YData',ydata);
    drawnow();

    %%calculate
    [IdIq] = calcTnPoint(motorModel,T_grid(ii),n_grid(ii));
    syreDriveSingt.Id = IdIq.Id;
    syreDriveSingt.Iq = IdIq.Iq;
    syreDriveSingt.T  = T_grid(ii);
    syreDriveSingt.n  = n_grid(ii);

    [i123] = MMM_eval_SyreDrivePoint(motorModel,syreDriveSingt);

    Iwave{ii}.time = linspace(0,60/n_grid(ii)/p,length(i123(1,:)));
    Iwave{ii}.Ia   = i123(1,:);
    Iwave{ii}.Ib   = i123(2,:);
    Iwave{ii}.Ic   = i123(3,:);
end

tic
cd(oldpath)

%% FEMM simulation - Iron/PM Loss
geo0 = motorModel.geo;
mat0 = motorModel.mat;
per0 = motorModel.per;

eval_type       = 'singmIron';
per0.custom_act = 1; % activate custom current .
per0.delta_sim_singt = 360/2;
per0.gamma = 0;
nPWM       = 10;     % 10 points each PWM step

Pfepm        = cell(Tstep,nstep);
Pfe_PWMfix   = zeros(Tstep,nstep);
Pfesh_PWMfix = zeros(Tstep,nstep);
Pfesc_PWMfix = zeros(Tstep,nstep);
Pferh_PWMfix = zeros(Tstep,nstep);
Pferc_PWMfix = zeros(Tstep,nstep);
Pfes_PWMfix  = zeros(Tstep,nstep);
Pfer_PWMfix  = zeros(Tstep,nstep);
Ppm_PWMfix   = zeros(Tstep,nstep);


parfor ii=1:numel(T_grid)
    %plot
    %     ii = vect(jj);
    xdata = [get(hfemm,'XData') n_grid(ii)];
    ydata = [get(hfemm,'YData') T_grid(ii)];
    set(hfemm,'XData',xdata,'YData',ydata);
    drawnow();

    %%calculate
    geoTmp = geo0;
    matTmp = mat0;
    perTmp = per0;
    perTmp.custom_ia  = Iwave{ii}.Ia;
    perTmp.custom_ib  = Iwave{ii}.Ib;
    perTmp.custom_ic  = Iwave{ii}.Ic;
    perTmp.EvalSpeed  = n_grid(ii);
    perTmp.nsim_singt = round((fPWM*nPWM)/(n_grid(ii)*p/60)/(360/perTmp.delta_sim_singt));
    if  perTmp.nsim_singt>1800
        perTmp.nsim_singt = 1800;
        disp('Check minimum speed: number of rotor positions exceeded the limit!')
    end
    [~,~,~,Pfepm{ii},~] = FEMMfitness([],geoTmp,perTmp,matTmp,eval_type,[pathname filename]);
    Pfe_PWMfix(ii)   = Pfepm{ii}.Pfe;
    Pfesh_PWMfix(ii) = Pfepm{ii}.Pfes_h;
    Pfesc_PWMfix(ii) = Pfepm{ii}.Pfes_c;
    Pferh_PWMfix(ii) = Pfepm{ii}.Pfer_h;
    Pferc_PWMfix(ii) = Pfepm{ii}.Pfer_c;
    Pfes_PWMfix(ii)  = Pfesh_PWMfix(ii)+Pfesc_PWMfix(ii);
    Pfer_PWMfix(ii)  = Pferh_PWMfix(ii)+Pferc_PWMfix(ii);
    Ppm_PWMfix(ii)   = Pfepm{ii}.Ppm;
end
t_FEMM = toc;
%% Copper Loss

R20 = Rs0/(1+0.004*(temp0-20));
Rs  = R20.*(1+0.004*(Tcu-20));

Pjs_PWMfix = zeros(Tstep,nstep);
Pjs_SIN = zeros(Tstep,nstep);

if fkac(1,end)<2*fPWM
    disp('AC Copper Loss: the max evaluated frequency is less than twice the PWM frequency! ')
    disp('Re-evaluate the kac for accurate results.')
end

parfor ii=1:numel(T_grid)
    Ia = Iwave{ii}.Ia;
    Fw = n_grid(ii)*p/60;

    L = length(Ia);
    Fs = Fw*L;
    I2 = abs(fft(Ia)/L);
    I1 = I2(1:floor(L/2+1));
    I1(2:end-1) = 2*I1(2:end-1);   %%CONTROLLARE IL X2
    fI = Fs*(0:(L/2))/L;

    fLoss = fI(fI<=max(max(fkac)));
    ILoss = I1(1:length(fLoss));

    kac_PWM = interp2(fkac,Tkac,kac,fLoss,Tcu.*ones(size(fLoss)));
    kac_PWM(isnan(kac_PWM)) = 1;
    Ploss = 3/2*Rs*(kac_PWM*l/(lend+l)+lend/(lend+l)).*ILoss.^2;
    Pjs_PWMfix(ii) = sum(Ploss);
    kac_SIN = interp2(fkac,Tkac,kac,Fw,Tcu);
    Pjs_SIN(ii) = 3/2*Rs*(kac_SIN*l/(lend+l)+lend/(lend+l))*max(ILoss)^2;
end
%% Sinusoidal Effy Map
motorModel.TnSetup.PMLossFactor   = 1;
motorModel.TnSetup.IronLossFactor = 1;
motorModel.TnSetup.IronLossFlag = 'Yes';
motorModel.TnSetup.PMLossFlag = 'Yes';

nVect = linspace(0,nmax*1.1,50);
TVect = linspace(0,Tmax*1.1,50);

%Initialize matrix results
[nmap,Tmap] = meshgrid(nVect,TVect);

TwMap.n     = nmap;             % speed reference [rpm]
TwMap.T     = Tmap;             % torque (mechanical) reference [Nm]
TwMap.Tem   = nan(size(nmap));  % electro-magnetic torque [Nm]
TwMap.Tout  = nan(size(nmap));  % torque produced (mechanical, for operating limits) [Nm]
TwMap.Id    = nan(size(nmap));  % d-axis magnetizing current [A]
TwMap.Iq    = nan(size(nmap));  % q-axis magnetizing current [A]
TwMap.Ploss = nan(size(nmap));  % Total loss [W]
TwMap.Pjs   = nan(size(nmap));  % Stator Joule loss (total) [W]
TwMap.PjDC  = nan(size(nmap));  % Stator Joule loss (only DC) [W]
TwMap.PjAC  = nan(size(nmap));  % Stator Joule loss (only AC) [W]
TwMap.Pfe   = nan(size(nmap));  % Total iron loss [W]
TwMap.Pfes  = nan(size(nmap));  % Stator iron loss [W]
TwMap.Pfer  = nan(size(nmap));  % Rotor iron loss [W]
TwMap.Ppm   = nan(size(nmap));  % Permanent magnet loss [W]
TwMap.Pjr   = nan(size(nmap));  % Rotor Joule loss [W]
TwMap.Pmech = nan(size(nmap));  % Mechanical loss [W]
TwMap.Ife   = nan(size(nmap));  % Current for iron loss [A]
TwMap.Rs    = nan(size(nmap));  % Stator resistance [Ohm]
TwMap.eff   = nan(size(nmap));  % Efficiency [pu]

for ii=1:numel(TwMap.n)
    [out] = calcTnPoint(motorModel,Tmap(ii),nmap(ii));
    if ~isnan(out.T)
        TwMap.Tout(ii)  = out.T;
        TwMap.Id(ii)    = out.Id;
        TwMap.Iq(ii)    = out.Iq;
        TwMap.Tem(ii)   = out.Tem;
        TwMap.P(ii)     = out.P;
        TwMap.Ploss(ii) = out.Ploss;
        TwMap.Pjs(ii)   = out.Pjs;
        TwMap.PjDC(ii)  = out.PjDC;
        TwMap.PjAC(ii)  = out.PjAC;
        TwMap.Pfe(ii)   = out.Pfe;
        TwMap.Pfes(ii)  = out.Pfes;
        TwMap.Pfer(ii)  = out.Pfer;
        TwMap.Ppm(ii)   = out.Ppm;
        TwMap.Pjr(ii)   = out.Pjr;
        TwMap.Pmech(ii) = out.Pmech;
        TwMap.Ife(ii)   = out.Ife;
        TwMap.Rs(ii)    = out.Rs;
        TwMap.eff(ii)   = out.eff;
    end
end
%% Interpolate and correction factors

Pfe_PWM  = griddata(n_grid,T_grid,Pfe_PWMfix,nmap,Tmap);
Pfes_PWM = griddata(n_grid,T_grid,Pfes_PWMfix,nmap,Tmap);
Pfer_PWM = griddata(n_grid,T_grid,Pfer_PWMfix,nmap,Tmap);
Ppm_PWM  = griddata(n_grid,T_grid,Ppm_PWMfix,nmap,Tmap);
Pjs_PWM  = griddata(n_grid,T_grid,Pjs_PWMfix,nmap,Tmap);

Pfe_PWM(isnan(TwMap.Pfe))  = nan;
Pfes_PWM(isnan(TwMap.Pfe)) = nan;
Pfer_PWM(isnan(TwMap.Pfe)) = nan;
Ppm_PWM (isnan(TwMap.Pfe)) = nan;
Pjs_PWM (isnan(TwMap.Pfe)) = nan;

kfe  = Pfe_PWM./TwMap.Pfe;
kfes = Pfes_PWM./TwMap.Pfes;
kfer = Pfer_PWM./TwMap.Pfer;
kpm  = Ppm_PWM./TwMap.Ppm;
kjs  = Pjs_PWM./TwMap.Pjs;

Ploss_PWM =  Pfe_PWM + Pjs_PWM + Ppm_PWM; %+ out.Pmech + out.Pjr;

TwMapPWM.n     = nmap;
TwMapPWM.Tout  = out.T;
TwMapPWM.T     = Tmap;
TwMapPWM.Id    = out.Id;
TwMapPWM.Iq    = out.Iq;
TwMapPWM.P     = TwMapPWM.n.*TwMapPWM.T*pi/30;
TwMapPWM.Ploss = Ploss_PWM;
TwMapPWM.Pjs   = Pjs_PWM;
TwMapPWM.PjDC  = out.PjDC;
TwMapPWM.PjAC  = Pjs_PWM-out.PjDC;
TwMapPWM.Pfe   = Pfe_PWM;
TwMapPWM.Pfes  = Pfes_PWM;
TwMapPWM.Pfer  = Pfer_PWM;
TwMapPWM.Ppm   = Ppm_PWM;
TwMapPWM.Pjr   = out.Pjr;
TwMapPWM.Pmech = out.Pmech;
TwMapPWM.Rs    = out.Rs;
TwMapPWM.eff   = TwMapPWM.P./(TwMapPWM.P+Ploss_PWM);

%% Output
resFolder = [motorModel.data.motorName '_results\MMM results\' 'TwMapPWM_' datestr(now,30) '\'];

kfix.kfe  = kfe;
kfix.kpm  = kpm;
kfix.kjs  = kjs;
kfix.kfes = kfes;
kfix.kfer = kfer;

TwMapSIN = TwMap;
TwMapPWM.n_grid = n_grid;
TwMapPWM.T_grid = T_grid;
TwMapSIN.nmap = nmap;
TwMapSIN.Tmap = Tmap;

if ~exist([pathname resFolder],'dir')
    mkdir([pathname resFolder]);
end
save([pathname resFolder 'TwMDataPWM.mat'],'motorModel','TwMapPWM','TwMapSIN','kfix');

MMM_MaxTwPWM_plot(kfix,motorModel,TwMapSIN,TwMapPWM,resFolder)   %%plot
