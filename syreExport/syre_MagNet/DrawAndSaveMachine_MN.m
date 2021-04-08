% Copyright 2018
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

function dataSet = DrawAndSaveMachine_MN(dataSet,filename,pathname)

% == NEW SCRIPT FOR SAVING A MACHINE FROM GUI MANUAL ENTRIES ==
% ==== TAKE THE DIR AND THE STRUCT DATASET ===============================

pathname = dataSet.currentpathname;
filename = dataSet.currentfilename;
load([pathname filename]);

eval_type = 'singt';

h = OpenMagnet(1);  % 1 = visible, 0 = invisible
draw_motor_in_MN(geo,mat,pathname,filename,h);
[h,f] = SaveDocumentMagnet(h,[pathname,filename(1:end-4),'.mn']);
CloseMagnet(h)

dataSet.slidingGap      = 1; % R347

% refresh GUI display data
load([pathname filename]);
dataSet.RQ = round(dataSet.RQ,4);
% dataSet.currentpathname = [pathname '\'];
dataSet.currentpathname = [pathname];
dataSet.currentfilename = filename;

