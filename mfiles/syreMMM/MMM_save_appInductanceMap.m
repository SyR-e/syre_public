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

function MMM_save_appInductanceMap(motorModel)

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'Apparent Inductance Maps - ' int2str(motorModel.data.tempPM) 'deg\'];

if ~exist([pathname resFolder],'dir')
    mkdir([pathname resFolder])
end

Id = motorModel.AppInductanceMap_dq.Id;
Iq = motorModel.AppInductanceMap_dq.Iq;
Ld = motorModel.AppInductanceMap_dq.Ld;
Lq = motorModel.AppInductanceMap_dq.Lq;
Fm = motorModel.AppInductanceMap_dq.Fm;


save([pathname resFolder 'appInductanceMap.mat'],'Id','Iq','Ld','Lq','Fm');