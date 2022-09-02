% Copyright 2020
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

function motorModel = MMM_CtrlSIM(motorModel)
    
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
    
    mainFolder = pwd;
    cd(ctrlFolder_path);         % Sets the correct "Current Folder" to run the simulation


    motorModel.SyreDrive.simOut = sim(motorModel.SyreDrive.SIM_path);   % Starts simulation
    cd(mainFolder);
    clear mex; close_system;
end
