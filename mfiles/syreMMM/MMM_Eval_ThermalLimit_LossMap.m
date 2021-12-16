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

function MMM_Eval_ThermalLimit_LossMap(motorModel,hax)

path      = motorModel.data.pathname;
file      = motorModel.data.motorName;
TCu_lim   = motorModel.Thermal.TempCuLimit;
TPm_lim   = motorModel.Thermal.TempPmLimit;
speed     = linspace(motorModel.Thermal.nmin,motorModel.Thermal.nmax,motorModel.Thermal.NumSpeed);

[filename, pathname] = uigetfile([path '\' file '_results\MMM results' '\*.mat'], 'Pick a Tw Map');
tmp = load([pathname filename]);

TwData = tmp.TwData;
TwMap  = tmp.TwMap;

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
plot(hax,TwMap.n(1,:),TwMap.T_top_W,'k','LineWidth',1.5);
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

interpWidth = 200;

T_m     = zeros(interpWidth,length(speed));
I_m     = zeros(interpWidth,length(speed));
Pfe_m   = zeros(interpWidth,length(speed));
Pfes_m  = zeros(interpWidth,length(speed));
Pfer_m  = zeros(interpWidth,length(speed));
Ppm_m   = zeros(interpWidth,length(speed));
Pjs_m   = zeros(interpWidth,length(speed));
Ploss_m = zeros(interpWidth,length(speed));

I = abs(TwMap.Id+j*TwMap.Iq);

for ii=1:length(speed)
    [~,index] = min(abs(TwMap.n(1,:)-speed(ii)));
    speed(ii) = TwMap.n(1,index);
    
    T_m(:,ii) = interp1(1:length(TwMap.T(2:end,index)),TwMap.T(2:end,index),linspace(1,length(TwMap.T(2:end,index)),interpWidth));
    I_m(:,ii) = interp1(1:length(I(2:end,index)),I(2:end,index),linspace(1,length(I(2:end,index)),interpWidth));
    
    Pfe_m(:,ii)   = interp1(1:length(TwMap.Pfe(2:end,index)),TwMap.Pfe(2:end,index),linspace(1,length(TwMap.Pfe(2:end,index)),interpWidth));
    Pfes_m(:,ii)  = interp1(1:length(TwMap.Pfes(2:end,index)),TwMap.Pfes(2:end,index),linspace(1,length(TwMap.Pfes(2:end,index)),interpWidth));
    Pfer_m(:,ii)  = interp1(1:length(TwMap.Pfer(2:end,index)),TwMap.Pfer(2:end,index),linspace(1,length(TwMap.Pfer(2:end,index)),interpWidth));
    Ppm_m(:,ii)   = interp1(1:length(TwMap.Ppm(2:end,index)),TwMap.Ppm(2:end,index),linspace(1,length(TwMap.Ppm(2:end,index)),interpWidth));
    Pjs_m(:,ii)   = interp1(1:length(TwMap.Pjs(2:end,index)),TwMap.Pjs(2:end,index),linspace(1,length(TwMap.Pjs(2:end,index)),interpWidth));
    Ploss_m(:,ii) = interp1(1:length(TwMap.Ploss(2:end,index)),TwMap.Ploss(2:end,index),linspace(1,length(TwMap.Ploss(2:end,index)),interpWidth));
   
    loss = [Pfe_m(:,ii) Pfes_m(:,ii) Pfer_m(:,ii) Ppm_m(:,ii) Pjs_m(:,ii) Ploss_m(:,ii)];
    
    [out] = MMM_Eval_ThermalLimit_LossMap_sing(mcad,loss,TCu_lim,TPm_lim,speed(ii));


    Ploss(ii) = out.Ploss;
    [~,index] = min(abs(Ploss_m(:,ii)-Ploss(ii)));
    if isfinite(Ploss(ii))
        T(ii)    = T_m(index,ii);
        I_th(ii) = I_m(index,ii);
    else
        T(ii)    = 0; 
        I_th(ii) = 0;
    end

    TempCu(ii)   = out.TempCu;
    TempPM(ii)   = out.TempPM;
    deltaWJ(ii) = out.deltaWJ;
    
    xdata(ii,:) = [speed(ii) speed(ii)];
    ydata(ii,:) = [0 T(ii)];
    
    plot(hax,xdata(ii,:),ydata(ii,:),'g','LineWidth',1.5);
    drawnow();
   
end

invoke(mcad,'Quit');

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

resFolder = [path file '_results\MMM results\Thermal\SteadyLossMap_Cu' num2str(TCu_lim) 'C_Pm' num2str(TPm_lim) 'C_' num2str(time)];
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
plot(TwMap.n(1,:),TwMap.T_top_W,'k');
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