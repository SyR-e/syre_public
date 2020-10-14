% Copyright 2018
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

function [SOL] = simulate_xdeg_MN(geo,per,eval_type,pathname,filename)

h = OpenMagnet(1);

switch eval_type
    case 'singt'
        nsim = per.nsim_singt;
        xdeg = per.delta_sim_singt;
        gamma = per.gamma;
    case 'singm'
        nsim = per.nsim_singt;
        xdeg = per.delta_sim_singt;
        gamma = per.gamma;
    otherwise
        error('MagNet Simulations not available for the selected evaluation type!');
end

% randFactor = 0;

th0 = geo.th0;
p   = geo.p;
r  = geo.r;
gap = geo.g;
% ns  = geo.ns;
pc  = 360/(6*geo.q*p)/2;
ps  = geo.ps;
n3ph = geo.win.n3phase; %AS number of 3-phase circuits
Q = 6*geo.q*geo.p*geo.win.n3phase;
skew_angle = 0;
N_parallel = 1;
N_cond = geo.win.Ns/geo.p/(geo.q)/size(geo.win.avv,1);
Br = per.BrPP;
tempPP = per.tempPP;
q = geo.q;
gamma_ = gamma*pi/180;

% evaluation of the phase current values for all positions to be simulated
iAmp = per.overload*per.i0;
% iAmp = per.overload*calc_io(geo,per);

id = iAmp * cos(gamma * pi/180);
iq = iAmp * sin(gamma * pi/180);
Imod = abs(id + 1i* iq);       % Modulo corrente (A)

% Avvio delle simulazioni e salvataggi
if  (eval_type == 'singt')
    
    iStr=num2str(Imod,3); iStr = strrep(iStr,'.','A');
    gammaStr=num2str(gamma,4); gammaStr = strrep(gammaStr,'.','d');
    if isempty(strfind(gammaStr, 'd'))
        gammaStr = [gammaStr 'd'];
    end
    nStr = int2str(per.EvalSpeed);
    nStr = strrep(nStr,'.','rpm');
    if ~strcmpi(nStr,'rpm')
        nStr = [nStr 'rpm'];
    end
    
    resFolder = [filename(1:end-4) '_results\FEA results\'];
    if ~exist([pathname resFolder],'dir')
        mkdir([pathname resFolder]);
    end
    
    CaseDir = [pathname resFolder filename(1:end-4) '_T_eval_',iStr,'_',gammaStr '_' nStr '_MN'];

else
    Idstr=num2str(id,3); Idstr = strrep(Idstr,'.','A');
    Iqstr=num2str(iq,3); Iqstr = strrep(Iqstr,'.','A');
    
    if isempty(strfind(Idstr, 'A'))
        Idstr = [Idstr 'A'];
    end
    if isempty(strfind(Iqstr, 'A'))
        Iqstr = [Iqstr 'A'];
    end
    
    CaseDir = [pathname 'CaseBackup\', filename(1:end-4),'_' Idstr 'x' Iqstr '_MN'];
    
end

mkdir(CaseDir);
CaseDir = [CaseDir '\'];

copyfile([pathname filename(1:end-4) '.mn'],[CaseDir filename(1:end-4) '.mn']);

DocSV = invoke(h.magnetHandler, 'openDocument',[CaseDir,filename(1:end-4) '.mn']);
Doc = invoke(h.magnetHandler, 'getDocument');
View = invoke(Doc, 'getCurrentView');
Solution = invoke(Doc, 'getSolution');
calculator = invoke (h.magnetHandler, 'getStackCalculator');

% set the operating parameters
if (Br==0)
    n_mag_simulati=0;
else
    %n_mag_simulati=size(geo.BLKLABELS.rotore.BarName,1);
    tmp=geo.BLKLABELS.rotore.xy(:,3);
    n_mag_simulati=length(tmp(tmp==6));
end

% type of solution (order)
% Polynomial Order = 1 (3 field points per triangle)
% Command = ['Call getDocument().setParameter("", "PolynomialOrder", "1", infoNumberParameter)'];
% invoke(MN6, 'processCommand', Command);

% SourcesOnAtTransientStart: smooth field at t = 0
Command = ['Call getDocument().setParameter("", "SourcesOnAtTransientStart", "Yes", infoStringParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

% 09 03 2010 - MeshAllowUnconstrainedHoles
Command = ['Call getDocument().setParameter("", "MeshAllowUnconstrainedHoles", "Yes", infoStringParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

% Model temperature
Command = ['Call getDocument().setParameter("", "Temperature", "' num2str(tempPP) '", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

% set transient options (time0 etc..)
time_end = 1000 * xdeg/(360*p)*60/per.EvalSpeed;   % ms
invoke(Doc,'setFixedIntervalTimeSteps',0, time_end, time_end/nsim);

% Iron loss time interval - global
Command = ['Call getDocument().setParameter("", "TransientAveragePowerLossStartTime", "0", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

Command = ['Call getDocument().setParameter("", "TransientAveragePowerLossEndTime", "' num2str(time_end) ' %ms", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

% Iron loss time interval - rotor
% time_start_pferot = time_end - time_end/Q;
%% temp -- macchina vagati
% time_start_pferot = time_end - q/Q*time_end;
time_start_pferot = 0;

Command = ['Call getDocument().setParameter("rotor", "TransientAveragePowerLossStartTime", "' num2str(time_start_pferot) ' %ms", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

Command = ['Call getDocument().setParameter("rotor", "TransientAveragePowerLossEndTime", "' num2str(time_end) ' %ms", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

Command = ['Call getDocument().setMotionSourceType("Motion#1", infoVelocityDriven)'];
invoke(h.magnetHandler, 'processCommand', Command);
%%

% Motion#1: transient angle Vs time
% start time = 0
Command = 'REDIM ArrayOfValues1(1)';
invoke(h.magnetHandler, 'processCommand', Command);
Command = 'ArrayOfValues1(0)= 0';
invoke(h.magnetHandler, 'processCommand', Command);
% end time == sim_angle (elt deg)
Command = ['ArrayOfValues1(1)= ' num2str(time_end)];
invoke(h.magnetHandler, 'processCommand', Command);

% start angle = Cas.skew_angle
Command = 'REDIM ArrayOfValues2(1)';
invoke(h.magnetHandler, 'processCommand', Command);
Command = ['ArrayOfValues2(0)= ' num2str(skew_angle)]; % Skew_New
invoke(h.magnetHandler, 'processCommand', Command);
% end angle = start + sim_angle
Command = ['ArrayOfValues2(1)= ' num2str(xdeg/p + skew_angle)];   %sim_angle elettrici % Skew_New
invoke(h.magnetHandler, 'processCommand', Command);
Command = 'Call getDocument().setMotionPositionVsTime("Motion#1", ArrayOfValues1, ArrayOfValues2)';
invoke(h.magnetHandler, 'processCommand', Command);

%%
phase_index = {'U','V','W'};
if exist('th0')
    % casi successivi feb 2010 (+90 perchè Waveform in Magnet\coil è un sin)
    phase_angle = [0 -120 120] + (gamma_*180/pi + th0(1) + 90);
else
    % casi precedenti feb 2010
    phase_angle = [0 -120 120] - 180/pi*(pi-gamma_);
end

for nn=1:n3ph
    for ii = 1:3
        Command=['Call getDocument().setCoilSourceType("',phase_index{ii},int2str(nn),'", infoCurrentDriven)'];
        invoke(h.magnetHandler, 'processCommand', Command);

        Command = 'REDIM ArrayOfValues(5)';
        invoke(h.magnetHandler, 'processCommand', Command);
        Command = 'ArrayOfValues(0)= 0';
        invoke(h.magnetHandler, 'processCommand', Command);
        Command = ['ArrayOfValues(1)= ' num2str(Imod/N_parallel)];
        invoke(h.magnetHandler, 'processCommand', Command);
        Command = ['ArrayOfValues(2)= ' num2str(per.EvalSpeed/60*p)];
        invoke(h.magnetHandler, 'processCommand', Command);
        Command = 'ArrayOfValues(3)= 0';
        invoke(h.magnetHandler, 'processCommand', Command);
        Command = 'ArrayOfValues(4)= 0';
        invoke(h.magnetHandler, 'processCommand', Command);
        Command = ['ArrayOfValues(5)= ' num2str(phase_angle(ii)+th0(nn)-th0(1))];
        invoke(h.magnetHandler, 'processCommand', Command);
        Command = ['Call getDocument.setSourceWaveform("' phase_index{ii} int2str(nn) '","SIN", ArrayOfValues)'];
        invoke(h.magnetHandler, 'processCommand', Command);
        % Command = ['getDocument.setParameter("' index{i} '", "WaveFormValues", "[0, ' num2str(Imod) ', ' num2str(n/60*p) ', 0, 0, 90]", infoArrayParameter)'];
        % invoke(h.magnetHandler, 'processCommand', Command);
        
        Command = ['Call getDocument().setCoilNumberOfTurns("U' int2str(nn) '", ' num2str(N_cond) ')'];
        invoke(h.magnetHandler, 'processCommand', Command);
        
        Command = ['Call getDocument().setCoilNumberOfTurns("V' int2str(nn) '", ' num2str(N_cond) ')'];
        invoke(h.magnetHandler, 'processCommand', Command);
        
        Command = ['Call getDocument().setCoilNumberOfTurns("W' int2str(nn) '", ' num2str(N_cond) ')'];
        invoke(h.magnetHandler, 'processCommand', Command);
        
    end
end

% %number of turns
% Command = ['Call getDocument().setCoilNumberOfTurns("U", ' num2str(N_cond) ')'];
% invoke(h.magnetHandler, 'processCommand', Command);
% Command = ['Call getDocument().setCoilNumberOfTurns("V", ' num2str(N_cond) ')'];
% invoke(h.magnetHandler, 'processCommand', Command);
% Command = ['Call getDocument().setCoilNumberOfTurns("W", ' num2str(N_cond) ')'];
% invoke(h.magnetHandler, 'processCommand', Command);

% solve
invoke(Doc, 'solveTransient2dWithMotion');

RunS_PostProcessingMN;

%save the model (.mn)
% Command=['Call getDocument().save("',CaseFileName,'")'];
Command=['Call getDocument().save("',[CaseDir,filename(1:end-4) '.mn'],'")'];
invoke(h.magnetHandler, 'processCommand', Command);

%% output
tempo   = time(2:end);         %ms
th = (geo.th0(1) * pi/180 + tempo/1000 * per.EvalSpeed * pi/30 * geo.p);

if (xdeg == 360)
    SOL.Pfes_h = Pfes_h;
    SOL.Pfes_c = Pfes_c;
    SOL.Pfer_h = Pfer_h;
    SOL.Pfer_c = Pfer_c;
    SOL.PjrBar = PjrBar;
    SOL.Ppm    = Ppm;
end
SOL.th = (th)';
SOL.id = tmp1(2:end)';
SOL.iq = tmp2(2:end)';
SOL.fd = Fluxd';
SOL.fq = Fluxq';
SOL.T = (mean([-Torque_1sim(:,1) Torque_1sim(:,2)],2))';

% Command='Call close(False)';
% invoke(h.magnetHandler, 'processCommand', Command);

CloseMagnet(h);

% cd(cd);






