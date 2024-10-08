% Copyright 2020
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

function [motorModelNew] = MMM_scale_allTemperatures(motorModelNew,motorModelOld)

motorModelNew = rmfield(motorModelNew,'PMtempModels');
motorModelNew.PMtempModels.tempVectPM = motorModelOld.PMtempModels.tempVectPM;
disp('Scaling maps at different PM temperature...')
for ii=1:length(motorModelNew.PMtempModels.tempVectPM)
    tmp = motorModelOld;
    tmp.FluxMap_dq       = motorModelOld.PMtempModels.FluxMap_dq{ii};
    tmp.FluxMap_dqt      = motorModelOld.PMtempModels.FluxMap_dqt{ii};
    tmp.IronPMLossMap_dq = motorModelOld.PMtempModels.IronPMLossMap_dq{ii};

    tmp = MMM_scale(tmp,motorModelNew.tmpScale);

    motorModelNew.PMtempModels.FluxMap_dq{ii}       = tmp.FluxMap_dq;
    motorModelNew.PMtempModels.FluxMap_dqt{ii}      = tmp.FluxMap_dqt;
    motorModelNew.PMtempModels.IronPMLossMap_dq{ii} = tmp.IronPMLossMap_dq;
    disp([' - map ' int2str(ii) ' of ' int2str(length(motorModelNew.PMtempModels.tempVectPM)) ' scaled'])
end
disp('All maps scaled!')


%% Old function
% pathnameNew  = motorModelNew.data.pathname;
% motorNameNew = motorModelNew.data.motorName;
% 
% pathnameOld  = motorModelOld.data.pathname;
% motorNameOld = motorModelOld.data.motorName;
% 
% tempVectPM = motorModelOld.data.tempVectPM;
% tempPMnew  = motorModelNew.data.tempPM;
% 
% clc
% disp('Scaling all the flux maps at different temperatures')
% 
% TempFolder = [pathnameNew motorNameNew '_results\MMM results\tempModels\'];
% if ~exist(TempFolder,'dir')
%     mkdir(TempFolder);
%     disp('Temperature models folder created')
% end
% 
% for ii=1:length(tempVectPM)
%     disp(['- Flux map ' int2str(ii) ' of ' int2str(length(tempVectPM)) ' - ' int2str(tempVectPM(ii)) 'deg'])
%     if tempVectPM(ii)==tempPMnew
%         disp([' flux map already scaled'])
%         motorModel = motorModelNew;
%     else
%         disp(['  scaling map...'])
%         load([pathnameOld motorNameOld '_results\MMM results\tempModels\motorModel_' int2str(tempVectPM(ii)) 'deg.mat'],'motorModel');
%         motorModel = MMM_scale(motorModel,motorModelNew.tmpScale);
%         motorModel.data.pathname  = pathnameNew;
%         motorModel.data.motorName = motorNameNew;
%         save([pathnameNew motorNameNew '_results\MMM results\tempModels\motorModel_' int2str(tempVectPM(ii)) 'deg.mat'],'motorModel')
%         disp(['  Flux map scaled and saved'])
%     end
% end
% 
% disp('Flux maps at different temperatures scaled and saved!')

