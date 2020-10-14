% Copyright 2020
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

function [motorModel] = MMM_defaultMotorModel()

data.pathname   = [fileparts(which('GUI_Syre.mlapp')) '\'];
data.motorName  = 'No motor selected';
data.i0         = 0;
data.Imax       = 1;
data.Vdc        = 0;
data.nmax       = 1;
data.n0         = 0;
data.T0         = 0;
data.P0         = 0;
data.p          = 2;
data.n3phase    = 1;
data.Rs         = 0;
data.tempCu     = 20;
data.tempPM     = [];
data.Ns         = 1;
data.l          = 0;
data.axisType   = 'SR';
data.motorType  = 'SR';
data.J          = 0;
data.tempVectPM = [];


data.nCurr = 4;

scaleFactors.Lld = 0;
scaleFactors.Llq = 0;
scaleFactors.Ns  = data.Ns;
scaleFactors.l   = data.l;

skewData.thSkw   = 0;
skewData.nSlice  = 1;
skewData.nPoints = 51;
% data.skew = skewData;

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
Tw.temperature      = 20;
Tw.MechLoss         = [0];
Tw.IronLossFlag     = 'No';
Tw.IronLossFactor   = 2;
Tw.SkinEffectFlag   = 'No';
Tw.SkinEffectMethod = 'LUT';
Tw.PMLossFlag       = 'No';
Tw.PMLossFactor     = 1;
Tw.Control          = 'Maximum efficiency';

SyreDrive.Ctrl_type = 'Current control';
SyreDrive.FMapsModel = 'dq Model';
SyreDrive.Converter.V0 = 0;
SyreDrive.Converter.Rd = 1e-4;
SyreDrive.Converter.dT = 1e-6;


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
motorModel.SyreDrive   = SyreDrive;
