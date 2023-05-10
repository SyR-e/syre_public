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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [motorModel] = MMM_updatePMtemperature(motorModel)

if isfield(motorModel,'PMtempModels')
    index = find(motorModel.PMtempModels.tempVectPM==motorModel.data.tempPM);
    if isempty(index)
        disp('PM temperature not present')
        motorModel.PMtempModels.tempVectPM = [motorModel.PMtempModels.tempVectPM motorModel.data.tempPM];
        index = length(motorModel.PMtempModels.tempVectPM);
        motorModel.PMtempModels.FluxMap_dq{index}       = motorModel.FluxMap_dq;
        motorModel.PMtempModels.FluxMap_dqt{index}      = motorModel.FluxMap_dqt;
        motorModel.PMtempModels.IronPMLossMap_dq{index} = motorModel.IronPMLossMap_dq;
        motorModel.data.tempVectPM = sort(motorModel.PMtempModels.tempVectPM);
    else
        motorModel.PMtempModels.FluxMap_dq{index}       = motorModel.FluxMap_dq;
        motorModel.PMtempModels.FluxMap_dqt{index}      = motorModel.FluxMap_dqt;
        motorModel.PMtempModels.IronPMLossMap_dq{index} = motorModel.IronPMLossMap_dq;
    end
else
    motorModel.PMtempModels.tempVectPM          = motorModel.data.tempPM;
    motorModel.PMtempModels.FluxMap_dq{1}       = motorModel.FluxMap_dq;
    motorModel.PMtempModels.FluxMap_dqt{1}      = motorModel.FluxMap_dqt;
    motorModel.PMtempModels.IronPMLossMap_dq{1} = motorModel.IronPMLossMap_dq;
end