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

function MMM_changeAxis_allTemperatures(motorModel)

pathname  = motorModel.data.pathname;
motorName = motorModel.data.motorName;

tempVectPM = motorModel.data.tempVectPM;

axisTypeNew = motorModel.data.axisType;

clc
disp('Changing the axis type at different temperatures')

TempFolder = [pathname motorName '_results\MMM results\tempModels\'];
if ~exist(TempFolder,'dir')
    mkdir(TempFolder);
    disp('Temperature models folder created')
end

for ii=1:length(tempVectPM)
    disp(['- Flux map ' int2str(ii) ' of ' int2str(length(tempVectPM)) ' - ' int2str(tempVectPM(ii)) 'deg'])
    disp(['  changing axis...'])
    load([pathname motorName '_results\MMM results\tempModels\motorModel_' int2str(tempVectPM(ii)) 'deg.mat'],'motorModel');
    motorModel = MMM_changeAxis(motorModel,axisTypeNew);
    motorModel.data.pathname  = pathname;
    motorModel.data.motorName = motorName;
    save([pathname motorName '_results\MMM results\tempModels\motorModel_' int2str(tempVectPM(ii)) 'deg.mat'],'motorModel')
    disp(['  Axis changed and model saved'])
end

disp('Models at different temperatures updated with new axis!')

