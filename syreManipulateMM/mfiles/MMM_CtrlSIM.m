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

    ctrlFolder_name = [motorModel.data.motorName '_ctrl'];
    ctrlFolder_path = [motorModel.data.pathname ctrlFolder_name];
    
    if ~exist(ctrlFolder_path,'dir')
        copyfile('D:\syre\Syr_ctrl_SIM', ctrlFolder_path);
    end
        
	save([ctrlFolder_path '\motorModel.mat'],'motorModel');
    
    MMM_print_MotorDataH(motorModel);
    % if MTPA / Flux Maps / Inductances have not been evaluated: evaluate them
    
    cd(ctrlFolder_path);                                                % Sets the correct "Current Folder
    
    motorModel.SyreDrive.SIM_path = [ctrlFolder_path '\Syr_ctrl_SFun.slx'];
    
    motorModel.SyreDrive.simOut = sim(motorModel.SyreDrive.SIM_path);   % Starts simulation
end
