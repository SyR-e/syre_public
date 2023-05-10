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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [motorModel,flag] = MMM_checkFile(pathname,filename)

if nargin()==0
    [filename,pathname] = uigetfile([cd '\*.mat'],'Select motor file');
    if ~filename
        error('No file selected')
    end
end


flag=0; % >1 if there are some errors

clc
disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
disp('Check motorModel structure for MMM')
disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
disp('Pathname:')
disp(['    ' pathname])
disp('Filename:')
disp(['    ' filename])

disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
disp('Chech motorModel')

load([pathname filename]);

if exist('motorModel','var')
    disp(' (v) motorModel loaded')
else
    disp(' (x) motorModel not present')
    error('Select a correct file!!!')
end


disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
disp('Chech pathname')
if strcmp(motorModel.data.pathname,pathname)
    disp(' (v) correct')
else
    motorModel.data.pathname = pathname;
    disp(' (x) wrong!!! Fixed')
    flag=1;
end

disp('Chech motorname')
if strcmp(motorModel.data.motorName,filename(1:end-4))
    disp(' (v) correct')
else
    motorModel.data.motorName = filename(1:end-4);
    disp(' (x) wrong!!! Fixed')
    flag=1;
end

disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')

disp('Check result folder')
resFolder = [pathname filename(1:end-4) '_results\'];
if exist(resFolder,'dir')
    disp(' (v) exist')
else
    mkdir(resFolder);
    disp(' (x) does not exist!!! Created')
    flag=1;
end

disp('Check MMM folder')
MMMfolder = [resFolder 'MMM results\'];
if exist(MMMfolder,'dir')
    disp(' (v) exist')
else
    mkdir(MMMfolder);
    disp(' (x) does not exist!!! Created')
    flag=1;
end

% disp('Check PM temperature models folder')
% tempFolder = [MMMfolder 'tempModels\'];
% if exist(tempFolder,'dir')
%     disp(' (v) exist')
% else
%     mkdir(tempFolder);
%     disp(' (x) does not exist!!! Created')
%     
% end

% disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
% 
% disp('Check flux maps at different temperatures')
% tempVectPM = motorModel.data.tempVectPM;
% disp([' Expected flux maps: ' int2str(length(tempVectPM))])
% unfound = 0;
% for ii=1:length(tempVectPM)
%     mapFile = [tempFolder 'motorModel_' int2str(tempVectPM(ii)) 'deg.mat'];
%     if exist(mapFile,'file')
%         disp([' (v) flux map at ' int2str(tempVectPM(ii)) 'deg found'])
%     else
%         disp([' (x) flux map at ' int2str(tempVectPM(ii)) 'deg not found'])
%         tempVectPM(ii) = NaN;
%         unfound = unfound+1;
%     end
% end
% 
% if unfound==0
%     disp(' (v) all the flux maps are present')
% else
%     disp([' (x) flux map cache corrupted. Missing ' int2str(unfound) ' maps!'])
%     motorModel.data.tempVectPM = tempVectPM(~isnan(tempVectPM));
%     flag=1;
% end
% 
% if sum(isnan(motorModel.data.tempVectPM))
%     motorModel.data.tempVectPM = motorModel.data.tempPM;
%     disp([' (v) added flux map at ' int2str(motorModel.data.tempPM) 'deg'])
% end
% 
% disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
% if prod(~isnan(motorModel.data.tempVectPM))
%     motorModel.data.tempVectPM = motorModel.data.tempPM;
% %     disp('Map added')
% end

if flag
    answer = 'Yes';
    answer = questdlg('Update file?','Update','Yes','No',answer);
    if strcmp(answer,'Yes')
        save([pathname filename],'motorModel','-append');
    end
end

