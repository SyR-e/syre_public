% Copyright 2020
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

function [motorModel] = MMM_back_compatibility(motorModel,Dflag)

if nargin()==1
    Dflag=1;
end

flag = 0;

if isfield(motorModel,'Tw')
    if ~isfield(motorModel.Tw,'PMLossFlag')
        motorModel.Tw.PMLossFlag = 'No';
        motorModel.Tw.PMLossFactor = 1;

        flag = 1;
        if Dflag
            disp('- Added PM loss and PM loss factor');
        end
    end
end

if ~isfield(motorModel.data,'J')
    if isfield(motorModel,'dataSet')
        if isfield(motorModel.dataSet,'RotorInertia')
            motorModel.data.J = motorModel.dataSet.RotorInertia;
        else
            motorModel.data.J = 0;
        end
    else
        motorModel.data.J = 0;
    end
    flag = 1;
    if Dflag
        disp('- Added rotor inertia');
    end
end

if ~isfield(motorModel,'SyreDrive')
    motorModel.SyreDrive.Ctrl_type = 'Torque control';
    motorModel.SyreDrive.FMapsModel = 'dq Model';
    motorModel.SyreDrive.Converter.V0 = 0;
    motorModel.SyreDrive.Converter.Rd = 1e-4;
    motorModel.SyreDrive.Converter.dT = 1;

    flag = 1;
    if Dflag
        disp('- Added Syre Drive');
    end
end

if ~isfield(motorModel.SyreDrive,'SS_on')
    motorModel.SyreDrive.SS_on = 'Off';
    motorModel.SyreDrive.SS_settings.inj_waveform = 'Sinusoidal';
    motorModel.SyreDrive.SS_settings.dem = 'Current';
    motorModel.SyreDrive.SS_settings.HS_ctrl = 'APP';

    flag = 1;
    if Dflag
        disp('- Added sensorless toggle');
    end
end

if ~isfield(motorModel.data,'tempVectPM')
    motorModel.data.tempVectPM = motorModel.data.tempPM;

    flag = 1;
    if Dflag
        disp('- Multiple PM temperature extension');
    end
end

if ~isfield(motorModel,'WaveformSetup')
    motorModel.WaveformSetup.CurrLoad  = motorModel.dqtElab.CurrLoad;
    motorModel.WaveformSetup.CurrAmpl  = motorModel.dqtElab.CurrAmpl;
    motorModel.WaveformSetup.CurrAngle = motorModel.dqtElab.CurrAngle;
    motorModel.WaveformSetup.EvalSpeed = motorModel.data.n0;
    motorModel.WaveformSetup.nCycle    = 1;
    %     motorModel = rmfield(motorModel,'dqtElab');
    flag = 1;
    if Dflag
        disp('- Added waveform tab and removed dqtMap tab');
    end
end

if ~isfield(motorModel.data,'R')
    if ~isempty(motorModel.dataSet)
        motorModel.data.R  = motorModel.dataSet.StatorOuterRadius;
        motorModel.scale.R = motorModel.dataSet.StatorOuterRadius;
    else
        motorModel.data.R = motorModel.data.l;
        motorModel.scale.R = motorModel.data.R;
    end
    flag = 1;
    if Dflag
        disp('- Added radial scaling');
    end
end

if ~isfield(motorModel.data,'lend')
    if ~isempty(motorModel.dataSet)
        if ~isfield(motorModel.dataSet,'EndWindingsLength')
            motorModel.dataSet.EndWindingsLength = calc_endTurnLength(motorModel.geo);
        end
        motorModel.data.lend  = motorModel.dataSet.EndWindingsLength;
    else
        motorModel.data.lend  = motorModel.data.l/3;
    end
    flag = 1;
    if Dflag
        disp('- Added end-winding input');
    end
end

% Update motorModel fields name and remove the wrong names
modFlag=0;

oldNames{1}  = 'fdfq';
oldNames{2}  = 'dqtMap';
oldNames{3}  = 'ironLoss';
oldNames{4}  = 'skinEffect';
oldNames{5}  = 'AOA';
oldNames{6}  = 'Inductance';
oldNames{7}  = 'idiq';
oldNames{8}  = 'dqtMapF';
oldNames{9}  = 'scale';
oldNames{10} = 'skew';
oldNames{11} = 'Tw';

newNames{1}  = 'FluxMap_dq';
newNames{2}  = 'FluxMap_dqt';
newNames{3}  = 'IronPMLossMap_dq';
newNames{4}  = 'acLossFactor';
newNames{5}  = 'controlTrajectories';
newNames{6}  = 'IncInductanceMap_dq';
newNames{7}  = 'FluxMapInv_dq';
newNames{8}  = 'FluxMapInv_dqt';
newNames{9}  = 'tmpScale';
newNames{10} = 'tmpSkew';
newNames{11} = 'TnSetup';

for ii=1:length(oldNames)
    if isfield(motorModel,oldNames{ii})
        command = ['motorModel.' newNames{ii} ' = motorModel.' oldNames{ii} ';'];
        eval(command);
        motorModel = rmfield(motorModel,oldNames{ii});
        modFlag=1;
    end
end

if modFlag
    if Dflag
        disp('- updated names of motorModel fields');
    end
    flag=1;
end


if ~isfield(motorModel,'Thermal')
    motorModel.Thermal.TempCuLimit = 180;
    motorModel.Thermal.TempPmLimit = 150;
    motorModel.Thermal.nmin        = 0;
    motorModel.Thermal.nmax        = motorModel.data.nmax;
    motorModel.Thermal.NumSpeed    = 6;
    flag = 1;
    if Dflag
        disp('- Added thermal section');
    end
end

% Custom geometry
if isfield(motorModel,'dataSet')
    if ~isfield(motorModel.dataSet,'custom')
        motorModel.dataSet.pShape.rotor  = [];
        motorModel.dataSet.pShape.stator = [];
        motorModel.dataSet.pShape.magnet = [];
        motorModel.dataSet.pShape.slot   = [];
        motorModel.dataSet.pShape.flag   = 0;
        motorModel.dataSet.custom        = 0;

        if Dflag
            disp('- Added Custom Geometry Mass Computation')
        end
        flag=1;
    end
end

% Thermal -Interpolate Flux Maps with PM temperature
if ~isfield(motorModel.Thermal,'interpTempPM')
    motorModel.Thermal.interpTempPM = 0;
    if Dflag
        disp('- Added Flux Maps Interpolation with PM temperature')
    end
    flag=1;
end

% Interpolate Flux Maps with PM temperature
if ~isfield(motorModel.data,'targetPMtemp')
    motorModel.data.targetPMtemp = 80;
    if Dflag
        disp('- Added Flux Maps Interpolation with PM temperature')
    end
    flag=1;
end

% Added model type in syreDrive (Andrei Bojoi MSc Thesis)
if ~isfield(motorModel.SyreDrive,'modelType')
    motorModel.SyreDrive.modelType = 'Average';
     if Dflag
        disp('- Added model type in syreDrive')
    end
    flag=1;
end

if ~isfield(motorModel.SyreDrive,'IronLoss')
    motorModel.SyreDrive.IronLoss = 'No';
    flag = 1;
    if Dflag
        disp('- Added IronLoss Syre Drive');
    end
end




% Improved PM temperature management (no cache, but all-in-one-file)
if ~isfield(motorModel,'PMtempModels')
    motorModel.PMtempModels.tempVectPM          = motorModel.data.tempPM;
    motorModel.PMtempModels.FluxMap_dq{1}       = motorModel.FluxMap_dq;
    motorModel.PMtempModels.FluxMap_dqt{1}      = motorModel.FluxMap_dqt;
    motorModel.PMtempModels.IronPMLossMap_dq{1} = motorModel.IronPMLossMap_dq;
    motorModel.data.tempVectPM = motorModel.data.tempPM;
    if Dflag
        disp('- Included PM temperature map in one file')
    end
    flag=1;
end

% Apparent inductance map
if ~isfield(motorModel,'AppInductanceMap_dq')
    motorModel.AppInductanceMap_dq = [];
    if Dflag
        disp('- Added apparent inductance map (dq)')
    end
    flag=1;
end

% PWM Effy Maps
if ~isfield(motorModel.SyreDrive.Converter,'fPWM')
    motorModel.TnSetup.Control = 'Max efficiency';
    motorModel.SyreDrive.Converter.fPWM = 10000;
    if Dflag
        disp('- Included PWM Efficiency Maps')
    end
    flag=1;
end

% Demagnetization
if ~isfield(motorModel,'DemagnetizationLimit')
    motorModel.DemagnetizationLimit = [];
    if Dflag
        disp('- Included Demagnetization Limit (current VS temperature)')
    end
    flag=1;
end

% Control strategy in syreDrive
if ~isfield(motorModel.SyreDrive,'Ctrl_strategy')
    motorModel.SyreDrive.Ctrl_strategy = 'FOC';
     if Dflag
        disp('- Included Control strategy in syreDrive')
    end
    flag=1;
end

if ~isfield(motorModel.TnSetup,'ASCsafeFlag')
   motorModel.TnSetup.ASCsafeFlag = 'No';
   motorModel.TnSetup.SkinEffectMethod = 'LUT';
   if Dflag
       disp('- Added the possibility to control with ASC safe state in torque-speed map computation')
   end
end


% message in command window if some data are added
if flag && Dflag
    msg = 'This project was created with an older version of SyR-e: proceed to SAVE MACHINE to update to latest version';
    %     title = 'WARNING';
    %     f = warndlg(msg,title,'modal');
    warning(msg);
end
