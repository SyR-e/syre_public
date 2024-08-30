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

function [motorModel] = MMM_load(pathname,filename)

mod = load([pathname filename]);

if isfield(mod,'motorModel')    % MMM already used 
    motorModel = mod.motorModel;
    motorModel.data.pathname = pathname;
    motorModel.data.motorName = filename(1:end-4);
elseif isfield(mod,'dataSet')   % first time MMM is used
    if isequal(mod.dataSet.currentfilename(1:end-4),filename(1:end-4))  % existing filename  
        data.motorName = mod.dataSet.currentfilename(1:end-4);
    else % new filename (e.g. a copy of another motor)
        data.motorName = filename(1:end-4);
        mod.dataSet.currentfilename = filename;
        dataSet = mod.dataSet;
        save([pathname filename],'dataSet','-append');
    end
    data.i0        = round(mod.dataSet.RatedCurrent,2);
    data.Imax      = round(mod.dataSet.RatedCurrent*mod.dataSet.CurrLoPP(1),2);
    data.Vdc       = 565;
    data.nmax      = mod.dataSet.OverSpeed;
    data.n0        = mod.dataSet.EvalSpeed;
    data.T0        = NaN;
    data.P0        = NaN;
    data.p         = mod.dataSet.NumOfPolePairs;
    data.n3phase   = mod.dataSet.Num3PhaseCircuit;
    data.Rs        = mod.dataSet.Rs;
    data.tempCu    = mod.dataSet.TargetCopperTemp;
    data.tempPM    = mod.dataSet.tempPP(1);
    data.Ns        = mod.dataSet.TurnsInSeries;
    data.l         = mod.dataSet.StackLength;
    data.lend      = mod.dataSet.EndWindingsLength;
    data.R         = mod.dataSet.StatorOuterRadius;
    if isfield(mod.dataSet,'RotorInertia')
        data.J = mod.dataSet.RotorInertia;
    else
        data.J = 0;
    end
    if(strcmp(mod.dataSet.TypeOfRotor,'Circular')||strcmp(mod.dataSet.TypeOfRotor,'Seg')||strcmp(mod.dataSet.TypeOfRotor,'ISeg')||strcmp(mod.dataSet.TypeOfRotor,'Fluid'))
        data.axisType = 'SR';
    elseif strcmp(mod.dataSet.TypeOfRotor,'IM')
        data.axisType = 'SR';
    else
        data.axisType = 'PM';
    end
    if strcmp(mod.dataSet.FluxBarrierMaterial,'Air')
        data.motorType = 'SR';
    else
        data.motorType = 'PM';
    end
    if strcmp(mod.dataSet.TypeOfRotor,'IM')
        data.motorType = 'IM';
    end
    
    data.pathname = pathname;
    
    data.Lld = 0;
    data.Llq = 0;
    
    data.nCurr = 1;
    
    data.tempVectPM = data.tempPM;
    
    skewData.thSkw   = 0;
    skewData.nSlice  = 1;
    skewData.nPoints = 51;
    
    scaleFactors.Lld = 0;
    scaleFactors.Llq = 0;
    scaleFactors.l   = data.l;
    scaleFactors.Ns  = data.Ns;
    scaleFactors.R   = data.R;
    
%     dqtElab.harmonic  = 6*[1 2 3];
%     dqtElab.CurrLoad  = 1;
%     dqtElab.CurrAmpl  = dqtElab.CurrLoad*data.i0;
%     dqtElab.CurrAngle = 45;
    
    Tw.nCurrent         = 1;
    Tw.nmin             = 0;
    Tw.nmax             = data.nmax;
    Tw.nstep            = 31;
    Tw.Tmin             = 0;
    Tw.Tmax             = 10;
    Tw.Tstep            = 31;
    Tw.temperature      = data.tempCu;
    Tw.MechLoss         = [0];
    Tw.IronLossFlag     = 'No';
    Tw.IronLossFactor   = 1;
    Tw.SkinEffectFlag   = 'No';
    Tw.SkinEffectMethod = 'LUT';
    Tw.PMLossFlag       = 'No';
    Tw.PMLossFactor     = 1;
    Tw.Control          = 'MTPA';
    
    SyreDrive.Ctrl_type    = 'Torque control';
    SyreDrive.FMapsModel   = 'dq Model';
    SyreDrive.Converter.V0 = 0;
    SyreDrive.Converter.Rd = 1e-4;
    SyreDrive.Converter.dT = 1;
    SyreDrive.SS_on = 'Off';
    SyreDrive.SS_settings.inj_waveform = 'Sinusoidal';
    SyreDrive.SS_settings.dem = 'Current';
    SyreDrive.SS_settings.HS_ctrl = 'APP';
    
    WaveformSetup.CurrLoad  = 1;
    WaveformSetup.CurrAmpl  = data.i0;
    WaveformSetup.CurrAngle = 45;
    WaveformSetup.EvalSpeed = data.n0;
    WaveformSetup.nCycle    = 1;
    
    motorModel.data                = data;
    motorModel.FluxMap_dq          = [];
    motorModel.FluxMap_dqt         = [];
    motorModel.IronPMLossMap_dq    = [];
    motorModel.acLossFactor        = [];
    motorModel.controlTrajectories = [];
    motorModel.IncInductanceMap_dq = [];
    motorModel.FluxMapInv_dq       = [];
    motorModel.FluxMapInv_dqt      = [];
    motorModel.tmpScale            = scaleFactors;
    motorModel.tmpSkew             = skewData;
    motorModel.TnSetup             = Tw;
    motorModel.SyreDrive           = SyreDrive;
    motorModel.WaveformSetup       = WaveformSetup;
    
    motorModel.dataSet = mod.dataSet;
    motorModel.geo     = mod.geo;
    motorModel.per     = mod.per;
    motorModel.mat     = mod.mat;
    
    motorModel.dataSet.currentpathname = pathname;
    motorModel.dataSet.currentfilename = [data.motorName '.mat'];
else
    error('Motor model not correct!')
end




