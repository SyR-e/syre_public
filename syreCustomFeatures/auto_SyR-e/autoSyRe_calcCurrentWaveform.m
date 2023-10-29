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

function autoSyRe_calcCurrentWaveform(dataSet)
%% Data
prompt = {'Reference torque [Nm]:','Reference speed [rpm]:','DC-link voltage [Vdc]:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'200','3000','300'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

T   = str2double(answer{1});
n   = str2double(answer{2});
Vdc = str2double(answer{3});

fPWM = 10000;   % To change the fPWM, it must be changed manually the User_Constants.h (#define fs 10.0e3f)

I     = nan;
gamma = nan;

% Vdc   = 400;
% T     = 150;
% n     = 1000;

%% Load motor
%[filename,pathname] = uigetfile('*.mat','Load machine model');
%load([dataSet.currentpathname dataSet.currentfilename]) %?

pathname = dataSet.currentpathname;
filename = dataSet.currentfilename;

motorModel = MMM_load(dataSet.currentpathname,dataSet.currentfilename);

%% Load flux-map

FEAfolder = [motorModel.data.pathname motorModel.data.motorName '_results\FEA results\'];
if ~exist(FEAfolder,'dir')
    FEAfolder = motorModel.data.pathname;
end

[filename1,pathname1] = uigetfile([FEAfolder '\*.mat'],'Load flux map');
data = load([pathname1 filename1]);
[motorModel.FluxMap_dq,~] = MMM_load_fdfq(data,motorModel.data.p);

%% AOA
method = 'LUT';
[AOA] = MMM_eval_AOA(motorModel,method);
motorModel.controlTrajectories = AOA;

%% Create Simulink Model
motorModel.SyreDrive.Converter.fPWM = fPWM;
motorModel.data.Vdc = Vdc;

motorModel.SyreDrive.modelType                = 'Istantaneous';
motorModel.SyreDrive.FMapsModel               = 'dq Model';
motorModel.SyreDrive.IronLoss                 = 'No';
% motorModel.SyreDrive.Ctrl_type                = 'Current control';
motorModel.SyreDrive.Ctrl_type                = 'Torque control';
% motorModel.SyreDrive.Ctrl_strategy            = 'FOC';
motorModel.SyreDrive.Ctrl_strategy            = 'DFVC';
motorModel.SyreDrive.SS_on                    = 'Off';
motorModel.SyreDrive.SS_settings.inj_waveform = 'Sinusoidal';
motorModel.SyreDrive.SS_settings.dem          = 'Current';
motorModel.SyreDrive.SS_settings.HS_ctrl      = 'Active Flux';

tic
if ~exist([pathname filename(1:(end-4)) '_ctrl_INST\'],'dir')
    motorModel = MMM_createSimulinkModel(motorModel);
else
    motorModel.SyreDrive.SIM_path = [pathname filename(1:(end-4)) '_ctrl_INST\' motorModel.data.motorName '_ctrl_INST.slx'];
end
time_buildSimulink = toc;

%% Run Simulink Model

if ~isnan(T)
    syreDriveSingt.Id = interp1(AOA.MTPA.T,AOA.MTPA.id,T);
    syreDriveSingt.Iq = interp1(AOA.MTPA.T,AOA.MTPA.iq,T);
else
    syreDriveSingt.Id = I*cosd(gamma);
    syreDriveSingt.Iq = I*sind(gamma);
end

syreDriveSingt.T    = T;
syreDriveSingt.n    = n;
syreDriveSingt.tSim = 0.45;

oldpath = pwd;
cd([motorModel.data.pathname motorModel.data.motorName '_ctrl_INST'])

tic
[i123] = MMM_eval_SyreDrivePoint(motorModel,syreDriveSingt);

time_SolveSimulink = toc;
cd(oldpath)

disp(['Build Simulink Model: ' num2str(time_buildSimulink) ' s'])
disp(['Solve Simulink Model: ' num2str(time_SolveSimulink) ' s'])

%% Save
Ia   = i123(1,:);
Ib   = i123(2,:);
Ic   = i123(3,:);
time = i123(4,:);

save([pathname motorModel.data.motorName '_results\FEA results\' '\customCurrent_' num2str(T) 'Nm_' num2str(n)  'rpm.mat'],'Ia','Ib','Ic','time');
disp(['Current waveform saved in: ' pathname motorModel.data.motorName '_results\FEA results\customCurrent_' num2str(T) 'Nm_' num2str(n)  'rpm.mat'])

