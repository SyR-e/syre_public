% Copyright 2023
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expres_trafs or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results] = skinEffect_evalCustomCurrent(dataSet)


%% FFT to select the number of steps
Ia   = dataSet.CustomCurrentA;
time = dataSet.CustomCurrentTime;

Fs = (time(2)-time(1))^-1;
L = length(Ia);
f_FFT = Fs*(0:(L/2))/L;
P2 = abs(fft(Ia)/L);
I_FFT = P2(1:floor(L/2+1))*2;

%fundamental
[~,ind] = max(I_FFT);
fmain = f_FFT(ind);
Imain = I_FFT(ind);


%% Load AC model
pathname = dataSet.currentpathname;
filename = dataSet.currentfilename;
FEAfolder = [pathname filename(1:(end-4)) '_results\FEA results\'];
if ~exist(FEAfolder,'dir')
    FEAfolder = pathname;
end

[filename1,pathname1] = uigetfile([FEAfolder '\*.mat'],'Load Skin Effect Results');
data = load([pathname1 filename1]);
% 
% if fkac(1,end)<2*fPWM
%     disp('AC Copper Loss: the max evaluated frequency is less than twice the PWM frequency! ')
%     disp('Re-evaluate the kac for accurate results.')
% end

fkac = data.results.f;
kac  = data.results.k;
Tkac = data.results.T;


%% Calculate
Tcu  = dataSet.TargetCopperTemp;
Rs   = dataSet.Rs;
l    = dataSet.StackLength;
lend = dataSet.EndWindingsLength;

if Tcu<min(min(Tkac)) || Tcu>max(max(Tkac)) 
    disp('AC Copper Loss: the target temperature is outside the evaluated temperatures boundaries! ')
end

kac_PWM = interp2(fkac,Tkac,kac,f_FFT,Tcu.*ones(size(f_FFT)));
kac_PWM(isnan(kac_PWM)) = 1;

Ploss = sum(3/2*Rs*(kac_PWM*l/(lend+l)+lend/(lend+l)).*I_FFT.^2);
PlossAC = 3/2*Rs.*(kac_PWM(ind)*l/(lend+l)+lend/(lend+l))*Imain.^2;
PlossDC = 3/2*Rs.*Imain.^2;

%% Output
mkdir([pathname1 'CustomCurrent'])

results.Ploss   = Ploss;
results.PlossAC = PlossAC;
results.PlossDC = PlossDC;


warning off

figure
figSetting(12,12,10)
subplot(2,1,1)
xlabel('$f$ [kHz]')
ylabel('$I$ [A]')
bar(f_FFT/1000,I_FFT,'b')
yyaxis right
plot(f_FFT/1000,kac_PWM)
xlim([0 max(max(fkac/1000))])
ylabel('$k_{AC}$')
subplot(2,1,2)
plot(time*1000,Ia)
xlabel('$t$ [ms]')
ylabel('$I$ [A]')
savefig([pathname1, 'CustomCurrent\Custom current.fig']);
savePrintFigure()

xNames{1} = '$PWM$';
xNames{2} = '$Fundamental$';
xNames{3} = '$DC$';

figure()
figSetting(12,12,10)
set(gca,'XLim',[0.5 3.5],'XTick',1:1:3,'XTickLabel',xNames);
b=bar([Ploss PlossAC PlossDC]);
xtips2 = b(1).XEndPoints;
ytips2 = b(1).YEndPoints;
labels2 = string(round(b(1).YData,2));
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
title('Copper Loss')
ylabel('$P_{Cu}$ [W]')
ax = gca;
ylim([0 ax.YTick(end)+ax.YTick(2)])

savefig([pathname1, 'CustomCurrent\Copper Loss.fig']);
savePrintFigure()

save([pathname1  'CustomCurrent\CopperLoss_Results.mat'],'Ploss','PlossAC','PlossDC','time');
disp(['Copper loss results saved in: ' pathname1 'CustomCurrent\'])