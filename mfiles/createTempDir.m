% Copyright 2014
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

function [thisfilepath,dirName]=createTempDir()

% createTempDir: creates a temporary variable for FEMM simulation

mang=num2str(randi((10^7),1));
dirName=mang(1:end);
warning off MATLAB:MKDIR:DirectoryExists
% thisfilepath = fileparts(which('MODEstart.m'));
thisfilepath = fileparts(which('GUI_Syre.mlapp'));
mkdir(fullfile(thisfilepath,'tmp'));
warning on MATLAB:MKDIR:DirectoryExists
while(exist([thisfilepath '\tmp\' dirName],'dir'))
    mang=num2str(randi((10^6),1));
    dirName=mang(1:end);
end
dirName=[thisfilepath '\tmp\' dirName '\'];
mkdir(dirName);
