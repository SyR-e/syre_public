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

function [motorModel] = MMM_changeAxis(motorModel,axisTypeNew)

if ~strcmp(axisTypeNew,motorModel.data.axisType)
    switch motorModel.data.axisType
        case 'SR'
            motorModel = MMM_sr2pm(motorModel);
        case 'PM'
            motorModel = MMM_pm2sr(motorModel);
    end
    
    % recompute AOA, Inductance, inverse models, if needed
    if ~isempty(motorModel.controlTrajectories)
        motorModel.controlTrajectories = MMM_eval_AOA(motorModel,motorModel.controlTrajectories.method);
    end
    
    if ~isempty(motorModel.IncInductanceMap_dq)
        motorModel.IncInductanceMap_dq = MMM_eval_inductanceMap(motorModel);
    end
    
    if ~isempty(motorModel.FluxMapInv_dq)
        motorModel.FluxMapInv_dq = MMM_eval_inverseModel_dq(motorModel);
    end
    
    if ~isempty(motorModel.FluxMapInv_dqt)
        motorModel.FluxMapInv_dqt = MMM_eval_inverse_dqtMap(motorModel);
    end
end

% check for other maps saved in the cache (other temperatures)

