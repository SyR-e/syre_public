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

function [motorModel] = MMM_changePMtemperature(motorModel,tempPM)

if strcmp(tempPM,'Add')
    % add new PM temperature
    FEAfolder = [motorModel.data.pathname motorModel.data.motorName '_results\FEA results\'];
    if ~exist('FEAfolder','dir')
        FEAfolder = motorModel.data.pathname;
    end
    [filename,pathname] = uigetfile([FEAfolder '\*.mat'],'Load new fdfq_idiq model');
    if filename
        data = load([pathname filename]);
        [fdfq,tempPM] = MMM_load_fdfq(data,motorModel.data.p);
        tempVect = motorModel.data.tempVectPM;
        if sum(tempVect==tempPM)
            msg = ['The maps at ' int2str(tempPM) ' deg are already loaded!'];
            f = warndlg(msg,'WARNING','modal');
        else
            indexOld = find(motorModel.PMtempModels.tempVectPM==motorModel.data.tempPM);
            motorModel.PMtempModels.FluxMap_dq{indexOld}       = motorModel.FluxMap_dq;
            motorModel.PMtempModels.FluxMap_dqt{indexOld}      = motorModel.FluxMap_dqt;
            motorModel.PMtempModels.IronPMLossMap_dq{indexOld} = motorModel.IronPMLossMap_dq;

            % loaded new PM temperature
            index = length(motorModel.PMtempModels.FluxMap_dq);
            motorModel.PMtempModels.tempVectPM(index+1) = tempPM;
            motorModel.PMtempModels.FluxMap_dq{index+1} = fdfq;
            % check for dqtMap
            if isfield(data,'dataSet')
                if data.dataSet.NumOfRotPosPP>=20
                    motorModel.PMtempModels.FluxMap_dqt{index+1} = MMM_eval_dqtMap(pathname,'F_map.mat');
                else
                    motorModel.PMtempModels.FluxMap_dqt{index+1} = [];
                end
            else
                motorModel.PMtempModels.FluxMap_dqt{index+1} = [];
            end
            % check for iron loss model
            motorModel.PMtempModels.IronPMLossMap_dq{index+1} = loadIronLossModel([pathname filename]);
            motorModel.data.tempVectPM = sort(motorModel.PMtempModels.tempVectPM);
            disp(['Added map at ' int2str(tempPM) 'deg'])


        end
    else
        error('Data not loaded')
    end
else
    % change PM temperature
    tempPM = eval(tempPM);
    
    indexOld = find(motorModel.PMtempModels.tempVectPM==motorModel.data.tempPM);
    motorModel.PMtempModels.FluxMap_dq{indexOld}       = motorModel.FluxMap_dq;
    motorModel.PMtempModels.FluxMap_dqt{indexOld}      = motorModel.FluxMap_dqt;
    motorModel.PMtempModels.IronPMLossMap_dq{indexOld} = motorModel.IronPMLossMap_dq;

    indexNew = find(motorModel.PMtempModels.tempVectPM==tempPM);
    motorModel.FluxMap_dq       = motorModel.PMtempModels.FluxMap_dq{indexNew};
    motorModel.FluxMap_dqt      = motorModel.PMtempModels.FluxMap_dqt{indexNew};
    motorModel.IronPMLossMap_dq = motorModel.PMtempModels.IronPMLossMap_dq{indexNew};
    motorModel.data.tempPM      = tempPM;
    if ~isempty(motorModel.controlTrajectories)
        motorModel.controlTrajectories = MMM_eval_AOA(motorModel,motorModel.controlTrajectories.method);
    end
    if ~isempty(motorModel.IncInductanceMap_dq)
        motorModel.IncInductanceMap_dq = MMM_eval_inductanceMap(motorModel);
    end
    if ~isempty(motorModel.AppInductanceMap_dq)
        motorModel.AppInductanceMap_dq = MMM_eval_appInductanceMap(motorModel);
    end
    if ~isempty(motorModel.FluxMapInv_dq)
        motorModel.FluxMapInv_dq = MMM_eval_inverseModel_dq(motorModel);
    end
    if ~isempty(motorModel.FluxMapInv_dqt)
        motorModel.FluxMapInv_dqt = MMM_eval_inverse_dqtMap(motorModel);
    end
end