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

function MMM_saveAs_allTemperatures(pathnameNew,filenameNew,pathnameOld,filenameOld)

if ~exist([pathnameOld filenameOld '_results\MMM results\tempModels\'],'dir')
    disp('tempModels not available with the original file')
else
    tempModelFiles = dir([pathnameOld filenameOld '_results\MMM results\tempModels\']);
    if length(tempModelFiles)>2
        disp('Saving tempModels files...')
        for ii=3:length(tempModelFiles)
            sourceFile = [pathnameOld filenameOld '_results\MMM results\tempModels\' tempModelFiles(ii).name];
            destinationFile = [pathnameNew filenameNew '_results\MMM results\tempModels\' tempModelFiles(ii).name];
            copyfile(sourceFile,destinationFile);
            disp(['- ' tempModelFiles(ii).name ' file copied'])
        end
        disp('tempModels saved!')
    else
        disp('No available tempModels')
    end
end