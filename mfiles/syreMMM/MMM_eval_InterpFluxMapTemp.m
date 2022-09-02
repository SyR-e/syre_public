% Copyright 2022
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

function [motorModel] = MMM_eval_InterpFluxMapTemp(motorModel)

path         = motorModel.data.pathname;
file         = motorModel.data.motorName;
targetPMtemp = motorModel.dataSet.targetPMtemp;

[filename1, pathname1] = uigetfile([path '\' file '_results\FEA results' '\*.mat'], 'Pick the low temperature Flux Map');
fdfq1 = load([pathname1 filename1]);
[filename2, pathname2] = uigetfile([path '\' file '_results\FEA results' '\*.mat'], 'Pick the high temperature Flux Map');
fdfq2 = load([pathname2 filename2]);

[fdfq] = interpFluxMapsTemperature(fdfq1,fdfq2,fdfq1.per.tempPP,fdfq2.per.tempPP,targetPMtemp);

motorModel.FluxMap_dq      = fdfq;
motorModel.data.tempVectPM = sort([motorModel.data.tempVectPM targetPMtemp]);
motorModel.data.tempPM     = targetPMtemp;

motorModel.FluxMap_dqt         = [];
motorModel.controlTrajectories = [];
motorModel.IncInductanceMap_dq = [];
motorModel.FluxMapInv_dq       = [];
motorModel.FluxMapInv_dqt      = [];

resFolder = [path file '_results\MMM results\tempModels\'];
mkdir(resFolder)
% motorModel = app.motorModel;
save([resFolder 'motorModel_' int2str(targetPMtemp) 'deg.mat'],'motorModel');