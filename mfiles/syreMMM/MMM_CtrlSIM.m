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

mainFolder = fileparts(which('GUI_Syre.mlapp'));
ctrlFolder_path = fileparts(motorModel.SyreDrive.SIM_path);

cd(ctrlFolder_path);         % Sets the correct "Current Folder" to run the simulation

motorModel.SyreDrive.simOut = sim(motorModel.SyreDrive.SIM_path);   % Starts simulation
cd(mainFolder);
clear mex; close_system;
end
