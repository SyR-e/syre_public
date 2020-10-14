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

function [motorModel] = MMM_back_compatibility(motorModel,Dflag)

if nargin()==1
    Dflag=1;
end

flag = 0;

if ~isfield(motorModel.Tw,'PMLossFlag')
    motorModel.Tw.PMLossFlag = 'No';
    motorModel.Tw.PMLossFactor = 1;
    
    flag = 1;
    if Dflag
        disp('- Added PM loss and PM loss factor');
    end
end

if ~isfield(motorModel.data,'J')
    if isfield(motorModel,'dataSet')
        if isfield(motorModel.dataSet,'RotorInertia')
            motorModel.data.J = motorModel.dataSet.RotorInertia;
        else
            motorModel.data.J = 0;
        end
    else
        motorModel.data.J = 0;
    end
    flag = 1;
    if Dflag
        disp('- Added rotor inertia');
    end
end

if ~isfield(motorModel,'SyreDrive')
    motorModel.SyreDrive.Ctrl_type = 'Current control';
    motorModel.SyreDrive.FMapsModel = 'dq Model';
    motorModel.SyreDrive.Converter.V0 = 0;
    motorModel.SyreDrive.Converter.Rd = 1e-4;
    motorModel.SyreDrive.Converter.dT = 1e-6;
    
    flag = 1;
    if Dflag
        disp('- Added Syre Drive');
    end
end

if ~isfield(motorModel.data,'tempVectPM')
    motorModel.data.tempVectPM = motorModel.data.tempPM;
    
    flag = 1;
    if Dflag
        disp('- Multiple PM temperature extension');
    end
end


% message in command window if some data are added
if flag && Dflag
    msg = 'This project was created with an older version of SyR-e: proceed to SAVE MACHINE to update to latest version';
    title = 'WARNING';
%     f = warndlg(msg,title,'modal');
    warning(msg);
end



