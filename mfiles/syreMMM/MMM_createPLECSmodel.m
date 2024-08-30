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



function motorModel = MMM_createPLECSmodel(motorModel)

% Check models
if isempty(motorModel.controlTrajectories)
    motorModel.controlTrajectories = MMM_eval_AOA(motorModel);
end

motorModel.FluxMapInv_dq = MMM_eval_inverse_dq_Simulink(motorModel); 

if isempty(motorModel.FluxMapInv_dqt)
    motorModel.FluxMapInv_dqt = MMM_eval_inverse_dqtMap(motorModel);
end


ctrlFolder_path = [motorModel.data.pathname motorModel.data.motorName '_ctrl_PLECS'];

syrePath = fileparts(which('GUI_Syre.mlapp'));

copyfile([syrePath '\syreDrive\PLECSmodel'], ctrlFolder_path);
movefile([ctrlFolder_path '\Motor_ctrl.plecs'],[ctrlFolder_path '\' motorModel.data.motorName '_Motor_ctrl.plecs']);

save([ctrlFolder_path '\motorModel.mat'],'motorModel');

MMM_print_MotorDataH_PLECS(motorModel);

motorModel.SyreDrive.SIM_path = [ctrlFolder_path '\' motorModel.data.motorName '_Model.plecs'];

disp('PLECS model created!')
disp(['pathname:'])
disp(['  ' ctrlFolder_path '\'])
disp(['filename:'])
disp(['  ' motorModel.data.motorName '_Motor_ctrl.plecs'])

run([ctrlFolder_path '\init_sim_PLECS.m'])




end