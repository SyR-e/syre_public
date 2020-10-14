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
    data.Imax      = round(mod.dataSet.RatedCurrent*mod.dataSet.CurrLoPP(1));
    data.Vdc       = 565;
    data.nmax      = mod.dataSet.OverSpeed;
    data.n0        = mod.dataSet.EvalSpeed;
    data.T0        = NaN;
    data.P0        = NaN;
    data.p         = mod.dataSet.NumOfPolePairs;
    data.n3phase   = mod.dataSet.Num3PhaseCircuit;
    data.Rs        = mod.dataSet.Rs;
    data.tempCu    = mod.dataSet.TargetCopperTemp;
    data.tempPM    = mod.dataSet.tempPP;
    data.Ns        = mod.dataSet.TurnsInSeries;
    data.l         = mod.dataSet.StackLength;
%     data.J         = mod.dataSet.RotorInertia;
    if(strcmp(mod.dataSet.TypeOfRotor,'Circular')||strcmp(mod.dataSet.TypeOfRotor,'Seg')||strcmp(mod.dataSet.TypeOfRotor,'ISeg')||strcmp(mod.dataSet.TypeOfRotor,'Fluid'))
        data.axisType = 'SR';
    else
        data.axisType = 'PM';
    end
    if strcmp(mod.dataSet.FluxBarrierMaterial,'Air')
        data.motorType = 'SR';
    else
        data.motorType = 'PM';
    end
    
    data.pathname = pathname;
    
    data.Lld = 0;
    data.Llq = 0;
    
    data.nCurr = 4;
    
    skewData.thSkw   = 0;
    skewData.nSlice  = 1;
    skewData.nPoints = 51;
    
    scaleFactors.Lld = 0;
    scaleFactors.Llq = 0;
    scaleFactors.l   = data.l;
    scaleFactors.Ns  = data.Ns;
    
    dqtElab.harmonic  = 6*[1 2 3];
    dqtElab.CurrLoad  = 1;
    dqtElab.CurrAmpl  = dqtElab.CurrLoad*data.i0;
    dqtElab.CurrAngle = 45;
    
    Tw.nCurrent         = 2;
    Tw.nmin             = 0;
    Tw.nmax             = data.nmax;
    Tw.nstep            = 11;
    Tw.Tmin             = 0;
    Tw.Tmax             = 10;
    Tw.Tstep            = 11;
    Tw.temperature      = data.tempCu;
    Tw.MechLoss         = [0];
    Tw.IronLossFlag     = 'No';
    Tw.IronLossFactor   = 2;
    Tw.SkinEffectFlag   = 'No';
    Tw.SkinEffectMethod = 'LUT';
    Tw.PMLossFlag       = 'No';
    Tw.PMLossFactor     = 1;
    Tw.Control          = 'Maximum efficiency';
    
    motorModel.data        = data;
    motorModel.fdfq        = [];
    motorModel.dqtMap      = [];
    motorModel.ironLoss    = [];
    motorModel.skinEffect  = [];
    motorModel.AOA         = [];
    motorModel.Inductance  = [];
    motorModel.idiq        = [];
    motorModel.dqtMapF     = [];
    motorModel.scale       = scaleFactors;
    motorModel.skew        = skewData;
    motorModel.dqtElab     = dqtElab;
    motorModel.Tw          = Tw;
    
    motorModel.dataSet = mod.dataSet;
    motorModel.geo     = mod.geo;
    motorModel.per     = mod.per;
    motorModel.mat     = mod.mat;
    
    motorModel.dataSet.currentpathname = pathname;
    motorModel.dataSet.currentfilename = [data.motorName '.mat'];
else
    error('Motor model not correct!')
end




