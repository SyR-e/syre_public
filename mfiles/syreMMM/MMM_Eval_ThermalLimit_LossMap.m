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

function MMM_eval_thermalLimit_LossMap(motorModel,hax)

path      = motorModel.data.pathname;
file      = motorModel.data.motorName;
TCu_lim   = motorModel.Thermal.TempCuLimit;
TPm_lim   = motorModel.Thermal.TempPmLimit;
speed     = linspace(motorModel.Thermal.nmin,motorModel.Thermal.nmax,motorModel.Thermal.NumSpeed);

if motorModel.Thermal.interpTempPM
    [filename1, pathname1] = uigetfile([path '\' file '_results\FEA results' '\*.mat'], 'Pick the low temperature Flux Map');
    fdfq1 = load([pathname1 filename1]);

    [filename2, pathname2] = uigetfile([path '\' file '_results\FEA results' '\*.mat'], 'Pick the high temperature Flux Map');
    fdfq2 = load([pathname2 filename2]);
end


[filename, pathname] = uigetfile([path '\' file '_results\MMM results' '\*.mat'], 'Pick a Tw Map');
tmp = load([pathname filename]);
TwData = tmp.TwData;
TwMap  = tmp.TwMap;
Tcu    = TwData.temperature;
clear tmp

% Check
if TwData.temperature ~= TCu_lim || motorModel.data.tempPM ~= TPm_lim
    disp('Warning: Torque map evaluated at different copper or magnet temperature')
end
if TwMap.Pfe(1,end)==0
    disp('Warning: The selected Torque map has no loss')
end

%% GUI axis
cla(hax)
set(hax,'XLim',[min(speed) max(speed)]);
set(hax,'YLim',[0 max(max(TwMap.T))]);
plot(hax,TwMap.n(1,:),TwMap.limits.Tmax,'k','LineWidth',1.5);
drawnow();

xdata = zeros(length(speed),2);
ydata = zeros(length(speed),2);

%% Initialization
T       = zeros(1,length(speed));
Ploss   = zeros(1,length(speed));
I_th    = zeros(1,length(speed));
TempCu  = zeros(1,length(speed));
TempPM  = zeros(1,length(speed));
deltaWJ = zeros(1,length(speed));

if nargin==1
    % no link to an existing axis...create a new axis
    figure()
    figSetting();
    hax = axes('OuterPosition',[0 0 1 1]);
end

%% Thermal simulations

filemot = [path file '.mot'];
mcad    = actxserver('MotorCAD.AppAutomation');
[fail]  = invoke(mcad,'LoadFromFile',filemot);
if fail == 1
    error('The requested MotorCAD model does not exist. Export the model first!')
end

interpWidthX = 1000;
interpWidthY = 1000;

T_m     = zeros(interpWidthX,interpWidthY);
I_m     = zeros(interpWidthX,interpWidthY);
Pfe_m   = zeros(interpWidthX,interpWidthY);
Pfes_m  = zeros(interpWidthX,interpWidthY);
Pfer_m  = zeros(interpWidthX,interpWidthY);
Ppm_m   = zeros(interpWidthX,interpWidthY);
Pjs_m   = zeros(interpWidthX,interpWidthY);
Ploss_m = zeros(interpWidthX,interpWidthY);

I = abs(TwMap.Id+1i*TwMap.Iq);
%%% interpolate maps
[X,Y]   = meshgrid(1:length(TwMap.T(1,:)),1:length(TwMap.T(:,1)));
[Xq,Yq] = meshgrid(linspace(1,length(X(1,:)),interpWidthX),linspace(1,length(Y),interpWidthY));

T_m     = interp2(X,Y,TwMap.T,Xq,Yq);
n       = interp2(X,Y,TwMap.n,Xq,Yq);
I_m     = interp2(X,Y,I,Xq,Yq);
Pfe_m   = interp2(X,Y,TwMap.Pfe,Xq,Yq);
Pfes_m  = interp2(X,Y,TwMap.Pfes,Xq,Yq);
Pfer_m  = interp2(X,Y,TwMap.Pfer,Xq,Yq);
Ppm_m   = interp2(X,Y,TwMap.Ppm,Xq,Yq);
Pjs_m   = interp2(X,Y,TwMap.Pjs,Xq,Yq);
Ploss_m = interp2(X,Y,TwMap.Ploss,Xq,Yq);


invoke(mcad,'SetVariable','StatorCopperLossesVaryWithTemp', 'True');
invoke(mcad,'SetVariable','StatorWindingTemperatureAtWhichPcuInput', Tcu);
%%%Setup the Magnetic Simulation to get the Magnet temperature (MCAD DOES NOT HAVE AN ACTIVEX PARAMETER FOR IT, hope they will add it)
% A magnetic simulation has to be executed with 1 points to optimize the
% computational time, the magnet temperature is passed to the magnetic
% module and then it can be read.
invoke(mcad,'SetVariable','TorqueCalculation','False');
invoke(mcad,'SetVariable','BackEMFCalculation','False');
invoke(mcad,'SetVariable','CoggingTorqueCalculation','False');
invoke(mcad,'SetVariable','TorqueSpeedCalculation','False');
invoke(mcad,'SetVariable','DemagnetizationCalc','False');
invoke(mcad,'SetVariable','ElectromagneticForcesCalc_Load','False');
invoke(mcad,'SetVariable','TorquePointsPerCycle',1);


for ii=1:length(speed)
    [~,index] = min(abs(n(1,:)-speed(ii)));
    speed(ii) = n(1,index);
    loss = [Pfe_m(:,index) Pfes_m(:,index) Pfer_m(:,index) Ppm_m(:,index) Pjs_m(:,index) Ploss_m(:,index)];

    [out] = MMM_eval_thermalLimit_LossMap_sing(mcad,loss,TCu_lim,TPm_lim,speed(ii));

    TempCu(ii)  = out.TempCu;
    TempPM(ii)  = out.TempPM;
    deltaWJ(ii) = out.deltaWJ;

    if motorModel.Thermal.interpTempPM
        [fdfq] = interpFluxMapsTemperature(fdfq1,fdfq2,fdfq1.per.tempPP,fdfq2.per.tempPP,TempPM(ii));
        tmp                     = motorModel;
        tmp.FluxMap_dq          = fdfq;
        tmp.TnSetup.temperature = TempCu(ii);

        Tref = linspace(0,max(TwMap.T(:,1)),100);

        Ploss_scaled = nan(1,length(Tref));
        T_scaled = nan(1,length(Tref));
        I_scaled = nan(1,length(Tref));

        for kk=1:length(Tref)
            [outTw] = calcTnPoint(tmp,Tref(kk),speed(ii));
            Ploss_scaled(kk) = outTw.Ploss;
            T_scaled(kk)     = outTw.T;
            I_scaled(kk)     = sqrt(outTw.Id.^2+outTw.Id.^2);
        end

        Ploss(ii) = out.Ploss;
        [~,index1] = min(abs(Ploss_scaled-Ploss(ii)));
        if isfinite(Ploss(ii))
            T(ii)    = T_scaled(index1);
            I_th(ii) = I_scaled(index1);
        else
            T(ii)    = 0;
            I_th(ii) = 0;
        end
    else
        Ploss(ii) = out.Ploss;
        [~,index1] = min(abs(loss(:,end)-Ploss(ii)));
        if isfinite(Ploss(ii))
            T(ii)    = T_m(index1,index);
            I_th(ii) = I_m(index1,index);
        else
            T(ii)    = 0;
            I_th(ii) = 0;
        end
    end

    xdata(ii,:) = [speed(ii) speed(ii)];
    ydata(ii,:) = [0 T(ii)];
    plot(hax,xdata(ii,:),ydata(ii,:),'g','LineWidth',1.5);
    drawnow();
end
% invoke(mcad,'Quit');

%% Output

P = T.*speed*pi/30;

ThTw.T       = T;
ThTw.P       = P;
ThTw.n       = speed;
ThTw.Loss    = Ploss;
ThTw.TempCu  = TempCu;
ThTw.TempPM  = TempPM;
ThTw.deltaWJ = deltaWJ;
ThTw.I       = I_th;

%% Plot
clock_1 = fix(clock);
time = join((string(clock_1(1:end-1))),"");

resFolder = [path file '_results\MMM results\Thermal\SteadyLossMap_Cu' num2str(TCu_lim) 'C_Pm' num2str(TPm_lim) 'C_' num2str(time) '\'];
mkdir(resFolder);

plot(hax,speed,T,'r','LineWidth',1.5);
drawnow();

hfig(1) = figure();
figSetting(15,10,12)
xlabel('$n$ [rpm]')
ylabel('$T$ [Nm]')
title('Continuos torque')
axis ([0 speed(end) 0 TwMap.T(end,1)])
set(hfig(1),'FileName',[resFolder 'TorqueVsSpeed.fig'])
plot(TwMap.n(1,:),TwMap.limits.Tmax,'k');
plot(speed,T,'b')
plot(speed,T,'ro','MarkerSize',4);

hfig(2) = figure();
figSetting(15,10,12)
xlabel('$n$ [rpm]')
ylabel('$P$ [W]')
title('Continuos power')
axis ([0 speed(end) 0 max(P)*1.2])
set(hfig(2),'FileName',[resFolder 'PowerVsSpeed.fig'])
plot(speed,P,'r')
plot(speed,P,'bo','MarkerSize',4)

%% Save figures

answer = 'No';
answer = questdlg('Save figures?','Save','Yes','No',answer);
if strcmp(answer,'Yes')
    if ~exist(resFolder,'dir')
        mkdir(resFolder);
    end

    for ii=1:length(hfig)
        savePrintFigure(hfig(ii));
    end

    save([resFolder 'ThTw_Data.mat'], 'ThTw');
end