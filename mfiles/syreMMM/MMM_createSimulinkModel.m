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

n_set = motorModel.data.n3phase;
motorModel.SyreDrive.Simulator = 'Simulink';

%% ------------------------- Compute Motor Maps----------------------------%
if isempty(motorModel.controlTrajectories)
    warning('Control trajectories not computed!')
end

motorModel.FluxMapInv_dq = MMM_eval_inverse_dq_Simulink(motorModel); 

if isempty(motorModel.FluxMapInv_dqt)
    motorModel.FluxMapInv_dqt = MMM_eval_inverse_dqtMap(motorModel);
end

if isempty(motorModel.IncInductanceMap_dq)
    motorModel.IncInductanceMap_dq = MMM_eval_inductanceMap(motorModel);
end

if isempty(motorModel.AppInductanceMap_dq)
    motorModel.AppInductanceMap_dq = MMM_eval_appInductanceMap(motorModel);
end

%% -----------------------------Simulink Generation-----------------------%

if(n_set>1)
    Generate_MultiThreePhase_Simulink(motorModel,n_set);
else

    % modelType = motorModel.SyreDrive.modelSetup.InverterModel;
    
    ctrlFolder_name = [motorModel.data.motorName '_ctrl'];
    
    ctrlFolder_path = [motorModel.data.pathname ctrlFolder_name];
    
    syrePath = fileparts(which('GUI_Syre.mlapp'));
    
            copyfile(checkPathSyntax([syrePath '\syreDrive\SimulinkModel']), ctrlFolder_path);
            movefile(checkPathSyntax([ctrlFolder_path '\Motor_ctrl.slx']),checkPathSyntax([ctrlFolder_path '\' motorModel.data.motorName '_ctrl.slx']));
    
    MMM_print_MotorDataH(motorModel);
    MMM_print_ConstantsH(motorModel);
    
    motorModel.SyreDrive.SIM_file = checkPathSyntax([ctrlFolder_path '\' motorModel.data.motorName '_ctrl.slx']);
    
    save(checkPathSyntax([ctrlFolder_path '\motorModel.mat']),'motorModel');
    
    disp('Simulink model created!')
    disp('pathname:')
    disp(['  ' ctrlFolder_path '\'])
    disp('filename:')
    disp(['  ' motorModel.data.motorName '_ctrl.slx'])


end


