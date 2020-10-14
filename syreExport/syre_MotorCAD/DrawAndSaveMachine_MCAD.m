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

function dataSet = DrawAndSaveMachine_MCAD(dataSet,filename,pathname)

% == NEW SCRIPT FOR SAVING A MACHINE FROM GUI MANUAL ENTRIES ==
% ==== TAKE THE DIR AND THE STRUCT DATASET ===============================

% dataSet = handles.dataSet;

filename=dataSet.currentfilename;
pathname=dataSet.currentpathname;
load([pathname filename]);

if strcmp(dataSet.TypeOfRotor,'Circular') || strcmp(dataSet.TypeOfRotor,'Seg') || strcmp(dataSet.TypeOfRotor,'SPM') || strcmp(dataSet.TypeOfRotor,'Vtype')
    draw_motor_in_MCAD(filename, pathname);
else
    error('Error: Export of the motor is not possible for the selected geometry.')
end

% dataSet.currentpathname = [pathname '\'];
dataSet.currentpathname = [pathname];
dataSet.currentfilename = filename;
end