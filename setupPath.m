% Copyright 2014
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

function setupPath(flagInfo)

if nargin==0
    flagInfo=0;
end

% add the required directories to the path
syreDirectory = fileparts(which('GUI_Syre.mlapp'));

% Handle platform-specific FEMM path (optional: skip or prompt user if not Windows)
if ispc
    addpath('C:\femm42\mfiles');
else
    warning('backtrace','off');
    warning('SyR-e is running on non-Windows system. FEMM not available!')
    warning('backtrace','on');
end

addpath(fullfile(syreDirectory));
addpath(fullfile(syreDirectory,'mfiles'));
addpath(genpath(fullfile(syreDirectory,'mfiles','MODE')));
addpath(genpath(fullfile(syreDirectory,'mfiles','syrmDesign')));
addpath(genpath(fullfile(syreDirectory,'mfiles','DemagAnalysis')));
addpath(genpath(fullfile(syreDirectory,'mfiles','MaterialLibraryFunctions')));
addpath(genpath(fullfile(syreDirectory,'mfiles','StructuralPDE')));
addpath(genpath(fullfile(syreDirectory,'mfiles','syreMMM')));
addpath(genpath(fullfile(syreDirectory,'mfiles','OctaveFunctions')));
addpath(genpath(fullfile(syreDirectory,'mfiles','ZenodoDatabaseFunctions')));

addpath(fullfile(syreDirectory,'materialLibrary'));
addpath(fullfile(syreDirectory,'motorExamples'));

addpath(genpath(fullfile(syreDirectory,'syreExport')));

% check additional features (custom functions)
customFeaturesDir = fullfile(syreDirectory,'syreCustomFeatures');
addpath(customFeaturesDir);
addon = dir(customFeaturesDir);
if length(addon)>2
    if flagInfo
        disp('Custom features added:')
    end
    for ii=3:length(addon)
        if addon(ii).isdir
            addpath(genpath(fullfile(customFeaturesDir, addon(ii).name)));
            if flagInfo
                disp(['- ' addon(ii).name]);
            end
        end
    end
end

% check for missing folders
resultsDir = fullfile(cd, 'results');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir)
end
tmpDir = fullfile(cd, 'tmp');
if ~exist(tmpDir,'dir')
    mkdir(tmpDir)
end

simMatFilesDir = fullfile(cd, 'syreDrive', 'PLECSModel', 'SimMatFiles');
if ~exist(simMatFilesDir,'dir')
    mkdir(simMatFilesDir)
end

% Check for custom library files
matLibDir = fullfile(syreDirectory, 'materialLibrary');
customFiles = {'custom_iron.mat', 'custom_layer.mat', 'custom_conductor.mat', 'custom_sleeve.mat'};
for ii = 1:length(customFiles)
    customFilePath = fullfile(matLibDir, customFiles{ii});
    if ~exist(customFilePath, 'file')
        MatLib = {};
        MatList = {};
        save(customFilePath, 'MatLib', 'MatList');
    end
end
