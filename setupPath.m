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

% clc

%if isoctave
%    thisfilepath = fileparts(canonicalize_file_name(thisfilepath));
%    endfileName = '\';
%end
addpath('C:\femm42\mfiles');
addpath(fullfile(syreDirectory));
addpath(fullfile(syreDirectory,'mfiles'));
addpath(fullfile(syreDirectory,'mfiles','MODE'));
addpath(fullfile(syreDirectory,'mfiles','syrmDesign'));
addpath(fullfile(syreDirectory,'mfiles','DemagAnalysis'));
addpath(fullfile(syreDirectory,'mfiles','MaterialLibraryFunctions'));
addpath(fullfile(syreDirectory,'mfiles','StructuralPDE'));
addpath(fullfile(syreDirectory,'mfiles','syreMMM'));
addpath(fullfile(syreDirectory,'mfiles','OctaveFunctions'));

addpath (fullfile(syreDirectory,'materialLibrary'));
addpath (fullfile(syreDirectory,'motorExamples'));

addpath (fullfile (syreDirectory,'syreExport'));
addpath(genpath(fullfile(syreDirectory,'syreExport\syre_Dxf')));
addpath(genpath(fullfile(syreDirectory,'syreExport\syre_MagNet')));
addpath(genpath(fullfile(syreDirectory,'syreExport\syre_MotorCAD')));
addpath(genpath(fullfile(syreDirectory,'syreExport\syre_AnsysMaxwell')));

% check additional features (custom functions)
addpath(fullfile(syreDirectory,'syreCustomFeatures'));
addon = dir([syreDirectory '\syreCustomFeatures\']);
if length(addon)>2
    if flagInfo
        disp('Custom features added:')
    end
    for ii=3:length(addon)
        addpath(genpath([syreDirectory '\syreCustomFeatures\' addon(ii).name]));
        if flagInfo
            disp(['- ' addon(ii).name]);
        end
    end
end

% savepath

% check for missing folders
if ~exist([cd '\results'],'dir')
    mkdir('results')
end
if ~exist([cd '\tmp'],'dir')
    mkdir('tmp')
end

if ~exist([cd '\syreDrive\PLECSModel\SimMatFiles'],'dir')
    mkdir('syreDrive\PLECSModel\SimMatFiles')
end
% Check for custom library files
if ~exist([syreDirectory '\materialLibrary\custom_iron.mat'],'file')
    MatLib = {};
    MatList = {};
    save([syreDirectory '\materialLibrary\custom_iron.mat'],'MatLib','MatList');
end
if ~exist([syreDirectory '\materialLibrary\custom_layer.mat'],'file')
    MatLib = {};
    MatList = {};
    save([syreDirectory '\materialLibrary\custom_layer.mat'],'MatLib','MatList');
end
if ~exist([syreDirectory '\materialLibrary\custom_conductor.mat'],'file')
    MatLib = {};
    MatList = {};
    save([syreDirectory '\materialLibrary\custom_conductor.mat'],'MatLib','MatList');
end
if ~exist([syreDirectory '\materialLibrary\custom_sleeve.mat'],'file')
    MatLib = {};
    MatList = {};
    save([syreDirectory '\materialLibrary\custom_sleeve.mat'],'MatLib','MatList');
end



