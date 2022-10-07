% Copyright 2022
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

function motorModel = MMM_createSimulinkModel(motorModel)

% Check models
if isempty(motorModel.controlTrajectories)
    motorModel.controlTrajectories = MMM_eval_AOA(motorModel);
end
if isempty(motorModel.FluxMapInv_dq)
    motorModel.FluxMapInv_dq = MMM_eval_inverseModel_dq(motorModel);
end
if isempty(motorModel.FluxMapInv_dqt)
    motorModel.FluxMapInv_dqt = MMM_eval_inverse_dqtMap(motorModel);
end
if isempty(motorModel.IncInductanceMap_dq)
    motorModel.IncInductanceMap_dq = MMM_eval_inductanceMap(motorModel);
end

modelType = motorModel.SyreDrive.modelType;

switch(modelType)
    case 'Average'
        ctrlFolder_name = [motorModel.data.motorName '_ctrl_AVG'];

    case 'Istantaneous'
        ctrlFolder_name = [motorModel.data.motorName '_ctrl_INST'];
end

ctrlFolder_path = [motorModel.data.pathname ctrlFolder_name];

syrePath = fileparts(which('GUI_Syre.mlapp'));

switch(modelType)
    case 'Average'
        copyfile([syrePath '\syreDrive\AVGModel'], ctrlFolder_path);
        movefile([ctrlFolder_path '\Motor_ctrl_AVG.slx'],[ctrlFolder_path '\' motorModel.data.motorName '_ctrl_AVG.slx']);
    case 'Istantaneous'
        copyfile([syrePath '\syreDrive\INSTModel'], ctrlFolder_path);
        movefile([ctrlFolder_path '\Motor_ctrl_INST.slx'],[ctrlFolder_path '\' motorModel.data.motorName '_ctrl_INST.slx']);
end

save([ctrlFolder_path '\motorModel.mat'],'motorModel');

MMM_print_MotorDataH(motorModel);

switch(modelType)
    case 'Average'
        motorModel.SyreDrive.SIM_path = [ctrlFolder_path '\' motorModel.data.motorName '_ctrl_AVG.slx'];

    case 'Istantaneous'
        motorModel.SyreDrive.SIM_path = [ctrlFolder_path '\' motorModel.data.motorName '_ctrl_INST.slx'];
end