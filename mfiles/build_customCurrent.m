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
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function dataSet = build_customCurrent(dataSet)

% pathname=dataSet.currentpathname;
% filemot = strrep(dataSet.currentfilename,'.mat','.fem');

if dataSet.EvalSpeed==0
    dataSet.EvalSpeed = 1000;
end
prompt = {'Phase Current [A]','Current phase angle [el deg]','Current ripple [A]','Rotor speed [rpm]','PWM frequency [kHz]'};
name   = 'Current Waveform Setup';
answer = {
    num2str(dataSet.RatedCurrent)
    num2str(dataSet.GammaPP)
    num2str(round(dataSet.RatedCurrent/20,2))
    num2str(dataSet.EvalSpeed)
    num2str(10)
    };
answer = inputdlg(prompt,name,1,answer);
Iamp   = eval(answer{1});
gamma  = eval(answer{2});
deltaI = eval(answer{3});
speed  = eval(answer{4});
fPWM   = eval(answer{5})*1000;

load([dataSet.currentpathname dataSet.currentfilename]);

th0 = geo.th0;
p = geo.p;

id = Iamp*cosd(gamma);
iq = Iamp*sind(gamma);
theta = linspace(th0,360+th0,1000)*pi/180;
time = linspace(0,60/speed/p,1000);
[iabc] = dq2abc(id,iq,theta);
ripple = 1/3*deltaI*cos(2*pi*fPWM*time)+2/3*deltaI*cos(2*pi*2*fPWM*time);
Ia = iabc(1,:)+ripple;
Ib = iabc(2,:)+ripple;
Ic = iabc(3,:)+ripple;

figure
figSetting(12,6,10)
plot(time*1000,Ia,'k')
plot(time*1000,Ib,'b')
plot(time*1000,Ic,'r')
xlabel('t [ms]')
ylabel('I [A]')
xlim([0 time(end)*1000])


%% Output
answer = 'No';
answer = questdlg('Save currents?','Save','Yes','No',answer);
if strcmp(answer,'Yes')
    [filename,pathname] = uiputfile(['newcustomcurrent.mat'],'input current name and location');
    save([pathname filename],'Ia','Ib','Ic','time');
    dataSet.CustomCurrentEnable = 1;
    dataSet.CustomCurrentA    = Ia;
    dataSet.CustomCurrentB    = Ib;
    dataSet.CustomCurrentC    = Ic;
    dataSet.CustomCurrentTime = time;
    dataSet.CustomCurrentFilename = filename;
end