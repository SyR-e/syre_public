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

function MMM_Eval_ThermalLimit_FluxMap(motorModel,hax)

path      = motorModel.data.pathname;
file      = motorModel.data.motorName;
TCu_lim   = motorModel.Thermal.TempCuLimit;
TPm_lim   = motorModel.Thermal.TempPmLimit;
I_th0     = motorModel.data.i0;
speed     = linspace(motorModel.Thermal.nmin,motorModel.Thermal.nmax,motorModel.Thermal.NumSpeed);

% Load Motor-CAD model
mcad        = actxserver('MotorCAD.AppAutomation');    %Open Motor-CAD
filenameMot = [file '.mot'];
[fail]      = invoke(mcad,'LoadFromFile',[path filenameMot]);
if fail == 1
    error('The requested MotorCAD model does not exist. Export the model first!')
end

%% Emag Simulation inputs

[filename, pathname] = uigetfile([path '\' file '_results\MMM results' '\*.mat'], 'Pick a Tw Map');
tmp=load([pathname filename]);
TwData = tmp.TwData;
TwMap  = tmp.TwMap;

I     = abs(TwMap.Id+j*TwMap.Iq);
Tem   = TwMap.Tem;
gamma = atand(TwMap.Iq./TwMap.Id);


if TwData.temperature ~= TCu_lim || motorModel.data.tempPM ~= TPm_lim
    disp('Warning: Torque map evaluated at different copper or magnet temperature')
end

tmp = TCu_lim;
invoke(mcad,'SetVariable','ArmatureConductor_Temperature',tmp);
tmp = TPm_lim;
invoke(mcad,'SetVariable','Magnet_Temperature',tmp);

%% Emag Simulation settings
invoke(mcad,'SetVariable','MagneticThermalCoupling',1);
invoke(mcad,'SetVariable','ProximityLossModel',0);
invoke(mcad,'SetVariable','MagneticSolver',1);

invoke(mcad,'SetVariable','TorqueCalculation','False');
invoke(mcad,'SetVariable','BackEMFCalculation','False');
invoke(mcad,'SetVariable','CoggingTorqueCalculation','False');
invoke(mcad,'SetVariable','TorqueSpeedCalculation','False');
invoke(mcad,'SetVariable','DemagnetizationCalc','False');
invoke(mcad,'SetVariable','ElectromagneticForcesCalc_Load','False');

invoke(mcad,'SetVariable','CurrentDefinition',0);


%% GUI axis
cla(hax)

set(hax,'XLim',[min(speed) max(speed)]);
set(hax,'YLim',[0 max(max(TwMap.T))]);
plot(hax,TwMap.n(1,:),TwMap.T_top_W,'k','LineWidth',1.5);
drawnow();

xdata = zeros(length(speed),2);
ydata = zeros(length(speed),2);

%% Initialization
T_m     = zeros(50,length(speed));
I_m     = zeros(50,length(speed));
gamma_m = zeros(50,length(speed));

T        = zeros(1,length(speed));
Loss_Tot = zeros(1,length(speed));
I_th     = zeros(1,length(speed));
TempCu   = zeros(1,length(speed));
TempPM   = zeros(1,length(speed));

%% Simulations

invoke(mcad,'ShowMagneticContext');

for ii=1:length(speed)
    [~,index] = min(abs(TwMap.n(1,:)-speed(ii)));
    speed(ii) = TwMap.n(1,index);
    tmp = speed(ii)+1;
    invoke(mcad,'SetVariable','ShaftSpeed',tmp);
    
    for kk=1:length(Tem(:,index))
        if isfinite(Tem(kk,index))
            T_m1(kk)     = Tem(kk,index);
            I_m1(kk)     = I(kk,index);
            gamma_m1(kk) = gamma(kk,index);
        else
            T_m     = interp1(1:length(T_m1(2:end)),T_m1(2:end),linspace(1,length(T_m1(2:end)),50));
            I_m     = interp1(1:length(I_m1(2:end)),I_m1(2:end),linspace(1,length(I_m1(2:end)),50));
            gamma_m = interp1(1:length(gamma_m1(2:end)),gamma_m1(2:end),linspace(1,length(gamma_m1(2:end)),50));
            break
        end
    end
    
    [~,jj] = min(abs(I_m-I_th0));
    I_sim  = I_m(jj);
    flag   = 1;
    flagPM = 0;
    
    while flag == 1
        flag      = 0;
        I_sim     = I_m(jj);
        gamma_sim = gamma_m(jj);
        invoke(mcad,'SetVariable','PeakCurrent',I_sim);
        invoke(mcad,'SetVariable','PhaseAdvance',gamma_sim);
        invoke(mcad,'SetVariable','MagneticThermalCoupling',1);
        invoke(mcad,'DoMagneticCalculation');
        invoke(mcad,'DoSteadyStateAnalysis');
        invoke(mcad,'SetVariable','MagneticThermalCoupling',2);
        invoke(mcad,'DoMagneticCalculation');
        invoke(mcad,'DoSteadyStateAnalysis');
        [~, WindingTemperature_Max] = invoke(mcad,'GetVariable','T_[Winding_Max]');
        [~, MagnetTemperature]      = invoke(mcad,'GetVariable','Magnet_Temperature');
        
        if WindingTemperature_Max<TCu_lim-10 && flagPM == 0
            jj = jj+1;
            flag  = 1;
        end
        if WindingTemperature_Max>TCu_lim+1 && flagPM == 0
            jj = jj-1;
            flag  = 1;
        end
        if MagnetTemperature>TPm_lim && flag == 0
            jj     = jj-1;
            flag   = 1;
            flagPM = 1;
        end
        
        if jj==1
            T(ii)      = nan;
            TempCu(ii) = nan;
            TempPM(ii) = nan;
            I_th(ii)   = nan;
            
            disp(['No continous point feasible at ' num2str(speed) ' rpm']);
            break
        end
    end
    
    if flag==0
        
        [~, WJ_out] = invoke(mcad,'GetVariable','WJ_Fluid_Outlet_Temp_[Active]');
        [~, WJ_in]  = invoke(mcad,'GetVariable','WJ_Fluid_Inlet_Temperature');
        deltaWJ = WJ_out - WJ_in;
        
        [~, Loss1] = invoke(mcad,'GetVariable','Power_Armature_Copper_Loss');
        [~, Loss2] = invoke(mcad,'GetVariable','Loss_[Stator_Back_Iron]');
        [~, Loss3] = invoke(mcad,'GetVariable','Loss_[Stator_Tooth]');
        [~, Loss4] = invoke(mcad,'GetVariable','Loss_[Magnet]');
        [~, Loss5] = invoke(mcad,'GetVariable','Loss_[Rotor_Back_Iron]');
        [~, Loss6] = invoke(mcad,'GetVariable','Loss_[Embedded_Magnet_Pole]');
        
        Loss_Tot(ii) = Loss1+Loss2+Loss3+Loss4+Loss5+Loss6;
        T(ii) = T_m(jj);
        TempCu(ii) = WindingTemperature_Max;
        TempPM(ii) = MagnetTemperature;
        I_th(ii)   = I_sim;
    end
    
    % plot GUI_MMM
    xdata(ii,:) = [speed(ii) speed(ii)];
    ydata(ii,:) = [0 T(ii)];
    plot(hax,xdata(ii,:),ydata(ii,:),'g','LineWidth',1.5);
    drawnow();
    
    % clear
    clear T_m1;   clear I_m1;   clear gamma_m1;
    clear T_m;    clear I_m;    clear gamma_m;
end

P = T.*speed*pi/30;

%% Output
out.T        = T;
out.P        = P;
out.n        = speed;
out.I        = I_th;
out.Loss     = Loss_Tot;
out.deltaTWJ = deltaWJ;
out.TempCu   = TempCu;
out.TempPM   = TempPM;

plot(hax,speed,T,'r','LineWidth',1.5);
drawnow();

%% Plot

clock_1 = fix(clock);
time = join((string(clock_1(1:end-1))),"");
resFolder = [path file '_results\MMM results\Thermal\SteadyFluxMap_Cu' num2str(TCu_lim) 'C_Pm' num2str(TPm_lim) 'C_' num2str(time)];
mkdir(resFolder);

hfig(1) = figure();
figSetting(11,6,8)
xlabel('$\omega$ [$rpm$]')
ylabel('$T$ [Nm]')
title('Continuos Torque')
axis ([0 speed(end) 0 max(T)*1.2])
set(hfig(1),'FileName',[resFolder '\TorqueVsSpeed.fig'])
plot(speed,T,'b')
plot(speed,T,'ro','MarkerSize',4)

hfig(2) = figure();
figSetting(11,6,8)
xlabel('$\omega$ [$rpm$]')
ylabel('$P$ [W]')
title('Continuos Power')
axis ([0 speed(end) 0 max(P)*1.2])
set(hfig(2),'FileName',[resFolder '\PowerVsSpeed.fig'])
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
    
    save([resFolder '\ThTw_Data.mat'], 'out');
end

