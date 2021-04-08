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

function [motorModel] = MMM_newModel()


[filename,pathname] = uiputfile([cd '\*.mat'],'Save new motorModel');


prompt = {...
    'Number of pole pairs',...          %  1
    'Number of three-phase sets',...    %  2
    'Motor type (SR/PM/IM)',...         %  3
    'Axis type (SR/PM)',...             %  4
    'Rated power [W]',...               %  5
    'Rated current [Apk]',...           %  6
    'DC link voltage [Vdc]',...         %  7
    'Rated speed [rpm]'...              %  8
    'Phase resistance [Ohm]',...        %  9
    'Winding temperature [°C]',...      % 10
    'Active length [mm]',...            % 11
    'Turns in series per phase',...     % 12
    'Skew angle',...                    % 13
    'Skew slices',...                   % 14
    'd-axis extra inductance [H]',...   % 15
    'q-axis extra inductance [H]',...   % 16
    };
name = 'Motor Ratings';
numlines = 1;
answers = {
    '2',...
    '1',...
    'SR',...
    'SR',...
    '1000',...
    '10',...
    '310',...
    '1500',...
    '0.2',...
    '20',...
    '100',...
    '120',...
    '0',...
    '1',...
    '0',...
    '0',...
    };
answers = inputdlg(prompt,name,numlines,answers);




data.pathname  = pathname;
data.motorName = filename(1:end-4);
data.i0        = eval(answers{6});
data.Imax      = 1.5*data.i0;
data.Vdc       = eval(answers{7});
data.n0        = eval(answers{8});
data.nmax      = 1.5*data.n0;
data.P0        = eval(answers{5});
data.T0        = round(data.P0/(data.n0*pi/30));
data.p         = eval(answers{1});
data.n3phase   = eval(answers{2});
data.Rs        = eval(answers{9});
data.tempCu    = eval(answers{10});
data.tempPM    = 20;
data.Ns        = eval(answers{12});
data.l         = eval(answers{11});
data.axisType  = answers{4};
data.motorType = answers{3};

data.nCurr = 4;

scaleFactors.Lld = eval(answers{15});
scaleFactors.Llq = eval(answers{16});
scaleFactors.Ns  = data.Ns;
scaleFactors.l   = data.l;

skewData.thSkw   = eval(answers{13});
skewData.nSlice  = eval(answers{14});
skewData.nPoints = 51;
% data.skew = skewData;

dqtElab.harmonic  = 6*[1 2 3];
dqtElab.CurrLoad  = 1;
dqtElab.CurrAmpl  = dqtElab.CurrLoad*data.i0;
dqtElab.CurrAngle = 45;

Tw.nmin             = 0;
Tw.nmax             = data.nmax;
Tw.nstep            = 11;
Tw.Tmin             = 0;
Tw.Tmax             = data.T0;
Tw.Tstep            = 11;
Tw.temperature      = data.tempCu;
Tw.MechLoss         = [0];
Tw.IronLossFlag     = 'No';
Tw.IronLossFactor   = 1.5;
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