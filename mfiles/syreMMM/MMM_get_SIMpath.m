% Copyright 2026
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

function [SIM_path] = MMM_get_SIMpath(motorModel)

pathname = motorModel.data.pathname;
motname  = motorModel.data.motorName;

if exist([pathname motname '_ctrl_INST\'],'dir')
    SIM_path =  [pathname motname '_ctrl_INST\' motname '_ctrl_INST.slx'];
    if ~exist(SIM_path,'file')
        SIM_path = [];
    end
else
    SIM_path = [];
end



