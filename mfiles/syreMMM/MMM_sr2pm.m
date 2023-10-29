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

function [motorModel] = MMM_sr2pm(motorModel)

% 1) dq maps
if ~isempty(motorModel.FluxMap_dq)
    fdfq = motorModel.FluxMap_dq;
    motorModel.FluxMap_dq.Id = fliplr(-fdfq.Iq');
    motorModel.FluxMap_dq.Iq = fliplr(+fdfq.Id');
    motorModel.FluxMap_dq.Fd = fliplr(-fdfq.Fq');
    motorModel.FluxMap_dq.Fq = fliplr(+fdfq.Fd');

    if isfield(fdfq,'T')
        motorModel.FluxMap_dq.T = fliplr(fdfq.T');
    end
    if isfield(fdfq,'dT')
        motorModel.FluxMap_dq.dT = fliplr(fdfq.dT');
    end
    if isfield(fdfq,'dTpp')
        motorModel.FluxMap_dq.dTpp = fliplr(fdfq.dTpp');
    end
end

% 2) dqtMaps
if ~isempty(motorModel.FluxMap_dqt)
    dqtMap = motorModel.FluxMap_dqt;
    motorModel.FluxMap_dqt.Id = sort(-dqtMap.Iq);
    motorModel.FluxMap_dqt.Iq = dqtMap.Id;
    
    motorModel.FluxMap_dqt.data.Id = flipud(pagetranspose(-dqtMap.data.Iq));
    motorModel.FluxMap_dqt.data.Iq = flipud(pagetranspose(dqtMap.data.Id));
    motorModel.FluxMap_dqt.data.Fd = flipud(pagetranspose(-dqtMap.data.Fq));
    motorModel.FluxMap_dqt.data.Fq = flipud(pagetranspose(dqtMap.data.Fd));
    %motorModel.FluxMap_dqt.data.T  = flipud(pagetranspose(dqtMap.data.T));
    %motorModel.FluxMap_dqt.data.th = flipud(pagetranspose(dqtMap.data.th));
    names = fieldnames(motorModel.FluxMap_dqt.data);
    for ii=1:length(names)
        if ~(strcmp(names{ii},'Id')||strcmp(names{ii},'Iq')||strcmp(names{ii},'Fd')||strcmp(names{ii},'Fq'))
            motorModel.FluxMap_dqt.data.(names{ii}) = flipud(pagetranspose(dqtMap.data.(names{ii})));
        end
    end
    
%     motorModel.FluxMap_dqt.fInt.Id = griddedInterpolant(motorModel.FluxMap_dqt.data.Id,motorModel.FluxMap_dqt.data.Iq,motorModel.FluxMap_dqt.data.th,motorModel.FluxMap_dqt.data.Id,'spline');
%     motorModel.FluxMap_dqt.fInt.Iq = griddedInterpolant(motorModel.FluxMap_dqt.data.Id,motorModel.FluxMap_dqt.data.Iq,motorModel.FluxMap_dqt.data.th,motorModel.FluxMap_dqt.data.Iq,'spline');
%     motorModel.FluxMap_dqt.fInt.Fd = griddedInterpolant(motorModel.FluxMap_dqt.data.Id,motorModel.FluxMap_dqt.data.Iq,motorModel.FluxMap_dqt.data.th,motorModel.FluxMap_dqt.data.Fd,'spline');
%     motorModel.FluxMap_dqt.fInt.Fq = griddedInterpolant(motorModel.FluxMap_dqt.data.Id,motorModel.FluxMap_dqt.data.Iq,motorModel.FluxMap_dqt.data.th,motorModel.FluxMap_dqt.data.Fq,'spline');
%     motorModel.FluxMap_dqt.fInt.T  = griddedInterpolant(motorModel.FluxMap_dqt.data.Id,motorModel.FluxMap_dqt.data.Iq,motorModel.FluxMap_dqt.data.th,motorModel.FluxMap_dqt.data.T,'spline');
end

% 3) ironLoss
if ~isempty(motorModel.IronPMLossMap_dq)
    ironLoss = motorModel.IronPMLossMap_dq;
    
    motorModel.IronPMLossMap_dq.Id     = fliplr(-ironLoss.Iq');
    motorModel.IronPMLossMap_dq.Iq     = fliplr(+ironLoss.Id');
    motorModel.IronPMLossMap_dq.Pfes_h = fliplr(ironLoss.Pfes_h');
    motorModel.IronPMLossMap_dq.Pfes_c = fliplr(ironLoss.Pfes_c');
    motorModel.IronPMLossMap_dq.Pfer_h = fliplr(ironLoss.Pfer_h');
    motorModel.IronPMLossMap_dq.Pfer_c = fliplr(ironLoss.Pfer_c');
    motorModel.IronPMLossMap_dq.Ppm    = fliplr(ironLoss.Ppm');
end

% update axisType
motorModel.data.axisType = 'PM';





