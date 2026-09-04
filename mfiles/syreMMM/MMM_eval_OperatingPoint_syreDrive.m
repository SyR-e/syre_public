% Copyright 2026
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

function [out] = MMM_eval_OperatingPoint_syreDrive(motorModel)

if ~isfield(motorModel.SyreDrive,'SIM_file')
    motorModel.SyreDrive.SIM_file = MMM_get_SIMpath(motorModel);
    if isempty(motorModel.SyreDrive.SIM_file)
        error('syreDrive model not available!')
    end
end

motorModel.SyreDrive.modelSetup.Ctrl_type     = 'Current control';
motorModel.SyreDrive.modelSetup.Ctrl_strategy = 'FOC';
motorModel.SyreDrive.modelSetup.InverterModel = 'PWM';

motorModel.FluxMapInv_dq = MMM_eval_inverse_dq_Simulink(motorModel);

% simulation
mainFolder = fileparts(which('GUI_Syre.mlapp'));
[ctrlFolder_path,modelName] = fileparts(motorModel.SyreDrive.SIM_file);
cd(ctrlFolder_path);         % Sets the correct "Current Folder" to run the simulation
save('motorModel.mat','motorModel')

Tpwm = 1/motorModel.SyreDrive.Converter.fPWM;

in = Simulink.SimulationInput(motorModel.SyreDrive.SIM_file);
in = setModelParameter(in,'StopTime','0.25');
% in = setModelParameter(in,'MaxStepSize',num2str(Tpwm/25));
in = in.setBlockParameter([modelName '/Load speed/LoadSpeedRamp'],'Slope','20');
in = in.setBlockParameter([modelName '/Load speed/LoadSpeedRamp'],'Start','0.05');

simOut = sim(in);

cd(mainFolder);
clear mex;
close_system(motorModel.SyreDrive.SIM_file,0);

% post-processing


t   = simOut.Outputs.id.Time;
T   = simOut.Out_M.T.Data;
n   = simOut.Out_M.n_m.Data;
id  = simOut.Outputs.id.Data;
iq  = simOut.Outputs.iq.Data;
fd  = simOut.Out_M.lambda_dq.Data(:,1);
fq  = simOut.Out_M.lambda_dq.Data(:,2);
iA  = simOut.Out_M.Iabc.Data(:,1);
iB  = simOut.Out_M.Iabc.Data(:,2);
iC  = simOut.Out_M.Iabc.Data(:,3);
thm = simOut.Out_M.theta_r.Data;

% sin_e = sin(thm*motorModel.data.p-motorModel.geo.th0*pi/180);
% cos_e = cos(thm*motorModel.data.p-motorModel.geo.th0*pi/180);
sin_e = sin(thm*motorModel.data.p);
cos_e = cos(thm*motorModel.data.p);
th = atan2(sin_e,cos_e)*180/pi;

fdq0 = [fd';fq';zeros(size(fd'))];

fabc = dq02abc(fdq0,th'*pi/180);

fA = fabc(1,:);
fB = fabc(2,:);
fC = fabc(3,:);

t1 = t(end)-60/n(end)*1/motorModel.data.p;

SOLraw.T  = T(t>=t1);
SOLraw.n  = n(t>=t1);
SOLraw.id = id(t>=t1);
SOLraw.iq = iq(t>=t1);
SOLraw.fd = fd(t>=t1);
SOLraw.fq = fq(t>=t1);
SOLraw.ia = iA(t>=t1);
SOLraw.ib = iB(t>=t1);
SOLraw.ic = iC(t>=t1);
SOLraw.fa = fA(t>=t1);
SOLraw.fb = fB(t>=t1);
SOLraw.fc = fC(t>=t1);
SOLraw.th = th(t>=t1);
SOLraw.t  = t(t>=t1);

SOLraw.th(SOLraw.th<0) = SOLraw.th(SOLraw.th<0)+360;

[SOLraw.th,index] = sort(SOLraw.th);
SOLraw.T  = SOLraw.T(index);
SOLraw.n  = SOLraw.n(index);
SOLraw.id = SOLraw.id(index);
SOLraw.iq = SOLraw.iq(index);
SOLraw.fd = SOLraw.fd(index);
SOLraw.fq = SOLraw.fq(index);
SOLraw.ia = SOLraw.ia(index);
SOLraw.ib = SOLraw.ib(index);
SOLraw.ic = SOLraw.ic(index);
SOLraw.fa = SOLraw.fa(index);
SOLraw.fb = SOLraw.fb(index);
SOLraw.fc = SOLraw.fc(index);
SOLraw.t  = SOLraw.t(index);

% re-interp over better grid for PWM
fPWM = motorModel.SyreDrive.Converter.fPWM;
p    = motorModel.geo.p;
n = SOLraw.n(end);

NUMSTEP = 2*fPWM/(n*p/60)*100;

SOL.T  = interp1(linspace(0,1,numel(SOLraw.T)),SOLraw.T,linspace(0,1,NUMSTEP));
SOL.n  = interp1(linspace(0,1,numel(SOLraw.T)),SOLraw.n,linspace(0,1,NUMSTEP));
SOL.id = interp1(linspace(0,1,numel(SOLraw.T)),SOLraw.id,linspace(0,1,NUMSTEP));
SOL.iq = interp1(linspace(0,1,numel(SOLraw.T)),SOLraw.iq,linspace(0,1,NUMSTEP));
SOL.fd = interp1(linspace(0,1,numel(SOLraw.T)),SOLraw.fd,linspace(0,1,NUMSTEP));
SOL.fq = interp1(linspace(0,1,numel(SOLraw.T)),SOLraw.fq,linspace(0,1,NUMSTEP));
SOL.ia = interp1(linspace(0,1,numel(SOLraw.T)),SOLraw.ia,linspace(0,1,NUMSTEP));
SOL.ib = interp1(linspace(0,1,numel(SOLraw.T)),SOLraw.ib,linspace(0,1,NUMSTEP));
SOL.ic = interp1(linspace(0,1,numel(SOLraw.T)),SOLraw.ic,linspace(0,1,NUMSTEP));
SOL.fa = interp1(linspace(0,1,numel(SOLraw.T)),SOLraw.fa,linspace(0,1,NUMSTEP));
SOL.fb = interp1(linspace(0,1,numel(SOLraw.T)),SOLraw.fb,linspace(0,1,NUMSTEP));
SOL.fc = interp1(linspace(0,1,numel(SOLraw.T)),SOLraw.fc,linspace(0,1,NUMSTEP));
SOL.th = interp1(linspace(0,1,numel(SOLraw.T)),SOLraw.th,linspace(0,1,NUMSTEP));
SOL.t  = interp1(linspace(0,1,numel(SOLraw.T)),SOLraw.t,linspace(0,1,NUMSTEP));


% out structure
out.id   = mean(SOL.id);
out.iq   = mean(SOL.iq);
out.fd   = mean(SOL.fd);
out.fq   = mean(SOL.fq);
out.T    = mean(SOL.T);
out.dT   = std(SOL.T);
out.dTpu = out.dT/out.T;
out.dTpp = max(abs(SOL.T))-min(abs(SOL.T));
out.IPF  = abs(sin(atan2(out.iq,out.id)-atan2(out.fq,out.fd)));
out.SOL  = SOL;

iAmp      = motorModel.WaveformSetup.CurrAmpl;
gamma     = motorModel.WaveformSetup.CurrAngle;
id = iAmp.*cosd(gamma);
iq = iAmp.*sind(gamma);
iStr=num2str(abs(id+j*iq),3);
iStr=strrep(iStr,'.','A');
gammaStr=num2str(atan2(iq,id)*180/pi,4);
gammaStr=strrep(gammaStr,'.','d');
if ~contains(gammaStr,'d')
    gammaStr=[gammaStr 'd'];
end


resFolder = checkPathSyntax([motorModel.data.motorName '_results\MMM results\' 'syreDrive_T_eval_' iStr '_' gammaStr ' - ' int2str(motorModel.data.tempPM) 'deg\']);
mkdir(motorModel.data.pathname,resFolder);

iabc = [SOL.ia;SOL.ib;SOL.ic;SOL.t;SOL.th];
Ia   = SOL.ia;
Ib   = SOL.ib;
Ic   = SOL.ic;
time = SOL.t;
th   = SOL.th;

% downsampling for FEA simulation.
nPoints = 30/(n*p)*10*2*fPWM*2;
% Downsample the data for FEA simulation
Ia   = interp1([1:1:numel(Ia  )]/numel(Ia  ),Ia  ,[1:1:nPoints]/nPoints);
Ib   = interp1([1:1:numel(Ib  )]/numel(Ib  ),Ib  ,[1:1:nPoints]/nPoints);
Ic   = interp1([1:1:numel(Ic  )]/numel(Ic  ),Ic  ,[1:1:nPoints]/nPoints);
time = interp1([1:1:numel(time)]/numel(time),time,[1:1:nPoints]/nPoints);
th   = interp1([1:1:numel(th  )]/numel(th  ),th  ,[1:1:nPoints]/nPoints);

save([motorModel.data.pathname resFolder 'out.mat'],'out','motorModel','iabc');
save([motorModel.data.pathname resFolder 'simOut.mat'],'simOut','-v7.3');
save([motorModel.data.pathname resFolder 'iabc.mat'],'Ia','Ib','Ic','time','th');

figure()
figSetting()
subplot(2,1,1)
plot(SOL.th,SOL.T)
title(['$T=' num2str(out.T) '\,Nm$']);
set(gca,'XLim',[0 360],'XTick',0:60:360)
xlabel('$\theta$ ($^\circ$)')
ylabel('(Nm)')
subplot(2,1,2)
plot(SOL.th,abs(sin(atan2(SOL.iq,SOL.id)-atan2(SOL.fq,SOL.fd))))
title(['$IPF=' num2str(out.IPF) '$']);
set(gca,'XLim',[0 360],'XTick',0:60:360)
xlabel('$\theta$ ($^\circ$)')
ylabel('IPF')

saveas(gcf,[motorModel.data.pathname resFolder 'singt_' iStr '_' gammaStr 'plot_Torque.fig']);

figure()
figSetting()
subplot(2,1,1)
plot(SOL.th,SOL.fd)
title(['$\lambda_d=' num2str(out.fd) 'Vs$'])
set(gca,'XLim',[0 360],'XTick',0:60:360)
xlabel('$\theta$ ($^\circ$)')
ylabel('(Vs)')
subplot(2,1,2)
plot(SOL.th,SOL.fq)
title(['$\lambda_q=' num2str(out.fq) 'Vs$'])
set(gca,'XLim',[0 360],'XTick',0:60:360)
xlabel('$\theta$ ($^\circ$)')
ylabel('(Vs)')

saveas(gcf,[motorModel.data.pathname resFolder 'singt_' iStr '_' gammaStr '_plot_Flux.fig'])

figure()
figSetting()
subplot(2,1,1)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ ($^\circ$)')
ylabel('$\lambda_{abc}$ (Vs)')
title(['Phase flux linkages'])
plot(SOL.th,SOL.fa);
plot(SOL.th,SOL.fb);
plot(SOL.th,SOL.fc);
subplot(2,1,2)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ ($^\circ$)')
ylabel('$i_{abc}$ (A)')
title(['Phase currents'])
plot(SOL.th,SOL.ia);
plot(SOL.th,SOL.ib);
plot(SOL.th,SOL.ic);

saveas(gcf,[motorModel.data.pathname resFolder 'singt_' iStr '_' gammaStr '_plot_Phase.fig'])
