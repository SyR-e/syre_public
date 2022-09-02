% Copyright 2021
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

function [motorModel] = MMM_pm2sr(motorModel)

% 1) dq maps
fdfq = motorModel.FluxMap_dq;
motorModel.FluxMap_dq.Id = flipud(fdfq.Iq');
motorModel.FluxMap_dq.Iq = flipud(-fdfq.Id');
motorModel.FluxMap_dq.Fd = flipud(fdfq.Fq');
motorModel.FluxMap_dq.Fq = flipud(-fdfq.Fd');

if isfield(fdfq,'T')
    motorModel.FluxMap_dq.T = flipud(fdfq.T');
end
if isfield(fdfq,'dT')
    motorModel.FluxMap_dq.dT = flipud(fdfq.dT');
end
if isfield(fdfq,'dTpp')
    motorModel.FluxMap_dq.dTpp = flipud(fdfq.dTpp');
end

% 2) dqtMaps
if ~isempty(motorModel.FluxMap_dqt)
    dqtMap = motorModel.FluxMap_dqt;
    motorModel.FluxMap_dqt.Id = dqtMap.Iq;
    motorModel.FluxMap_dqt.Iq = sort(-dqtMap.Id);

    motorModel.FluxMap_dqt.data.Id = fliplr(pagetranspose(dqtMap.data.Iq));
    motorModel.FluxMap_dqt.data.Iq = fliplr(-pagetranspose(dqtMap.data.Id));
    motorModel.FluxMap_dqt.data.Fd = fliplr(pagetranspose(dqtMap.data.Fq));
    motorModel.FluxMap_dqt.data.Fq = fliplr(-pagetranspose(dqtMap.data.Fd));
    motorModel.FluxMap_dqt.data.T  = fliplr(pagetranspose(dqtMap.data.T));
    motorModel.FluxMap_dqt.data.th = fliplr(pagetranspose(dqtMap.data.th));

    motorModel.FluxMap_dqt.fInt.Id = griddedInterpolant(motorModel.FluxMap_dqt.data.Id,motorModel.FluxMap_dqt.data.Iq,motorModel.FluxMap_dqt.data.th,motorModel.FluxMap_dqt.data.Id,'spline');
    motorModel.FluxMap_dqt.fInt.Iq = griddedInterpolant(motorModel.FluxMap_dqt.data.Id,motorModel.FluxMap_dqt.data.Iq,motorModel.FluxMap_dqt.data.th,motorModel.FluxMap_dqt.data.Iq,'spline');
    motorModel.FluxMap_dqt.fInt.Fd = griddedInterpolant(motorModel.FluxMap_dqt.data.Id,motorModel.FluxMap_dqt.data.Iq,motorModel.FluxMap_dqt.data.th,motorModel.FluxMap_dqt.data.Fd,'spline');
    motorModel.FluxMap_dqt.fInt.Fq = griddedInterpolant(motorModel.FluxMap_dqt.data.Id,motorModel.FluxMap_dqt.data.Iq,motorModel.FluxMap_dqt.data.th,motorModel.FluxMap_dqt.data.Fq,'spline');
    motorModel.FluxMap_dqt.fInt.T  = griddedInterpolant(motorModel.FluxMap_dqt.data.Id,motorModel.FluxMap_dqt.data.Iq,motorModel.FluxMap_dqt.data.th,motorModel.FluxMap_dqt.data.T,'spline');
end

% 3) ironLoss
if ~isempty(motorModel.IronPMLossMap_dq)
    ironLoss = motorModel.IronPMLossMap_dq;

    motorModel.IronPMLossMap_dq.Id     = flipud(ironLoss.Iq');
    motorModel.IronPMLossMap_dq.Iq     = flipud(-ironLoss.Id');
    motorModel.IronPMLossMap_dq.Pfes_h = flipud(ironLoss.Pfes_h');
    motorModel.IronPMLossMap_dq.Pfes_c = flipud(ironLoss.Pfes_c');
    motorModel.IronPMLossMap_dq.Pfer_h = flipud(ironLoss.Pfer_h');
    motorModel.IronPMLossMap_dq.Pfer_c = flipud(ironLoss.Pfer_c');
    motorModel.IronPMLossMap_dq.Ppm    = flipud(ironLoss.Ppm');
end

% update axisType
motorModel.data.axisType = 'SR';





