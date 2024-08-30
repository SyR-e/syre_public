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

%% ------------------------- Compute Motor Maps----------------------------%
if isempty(motorModel.controlTrajectories)
    % motorModel.controlTrajectories = MMM_eval_AOA(motorModel,'LUT');
    warning('Control trajectories not computed!')
end

motorModel.FluxMapInv_dq = MMM_eval_inverse_dq_Simulink(motorModel); 

if isempty(motorModel.FluxMapInv_dqt)
    % motorModel.FluxMapInv_dqt = MMM_eval_inverse_dqtMap(motorModel);
    warning('dqt inverse maps model not computed!')
end

if isempty(motorModel.IncInductanceMap_dq)
    % motorModel.IncInductanceMap_dq = MMM_eval_inductanceMap(motorModel);
    warning('Incremental inductance maps not computed!')
end

%% -----------------------------Simulink Generation-----------------------%

if(n_set>1)
    Generate_MultiThreePhase_Simulink(motorModel,n_set);
else

    modelType = motorModel.SyreDrive.modelSetup.modelType;
    
    switch(modelType)
        case 'Average'
            % ctrlFolder_name = [motorModel.data.motorName '_ctrl_AVG'];
            ctrlFolder_name = [motorModel.data.motorName '_ctrl_INST'];
    
        case 'Istantaneous'
            ctrlFolder_name = [motorModel.data.motorName '_ctrl_INST'];
    end
    
    ctrlFolder_path = [motorModel.data.pathname ctrlFolder_name];
    
    syrePath = fileparts(which('GUI_Syre.mlapp'));
    
    switch(modelType)
        case 'Average'
            % copyfile([syrePath '\syreDrive\AVGModel'], ctrlFolder_path);
            % movefile([ctrlFolder_path '\Motor_ctrl_AVG.slx'],[ctrlFolder_path '\' motorModel.data.motorName '_ctrl_AVG.slx']);
            copyfile([syrePath '\syreDrive\INSTModel'], ctrlFolder_path);
            movefile([ctrlFolder_path '\Motor_ctrl_INST.slx'],[ctrlFolder_path '\' motorModel.data.motorName '_ctrl_INST.slx']);
        case 'Istantaneous'
            copyfile([syrePath '\syreDrive\INSTModel'], ctrlFolder_path);
            movefile([ctrlFolder_path '\Motor_ctrl_INST.slx'],[ctrlFolder_path '\' motorModel.data.motorName '_ctrl_INST.slx']);
    end
    
    
    MMM_print_MotorDataH(motorModel);
    
    switch(modelType)
        case 'Average'
            % motorModel.SyreDrive.SIM_path = [ctrlFolder_path '\' motorModel.data.motorName '_ctrl_AVG.slx'];
            motorModel.SyreDrive.SIM_path = [ctrlFolder_path '\' motorModel.data.motorName '_ctrl_INST.slx'];
        case 'Istantaneous'
            motorModel.SyreDrive.SIM_path = [ctrlFolder_path '\' motorModel.data.motorName '_ctrl_INST.slx'];
    end
    
    save([ctrlFolder_path '\motorModel.mat'],'motorModel');
    
    disp('Simulink model created!')
    disp(['pathname:'])
    disp(['  ' ctrlFolder_path '\'])
    disp(['filename:'])
    disp(['  ' motorModel.data.motorName '_ctrl_INST.slx'])


end


