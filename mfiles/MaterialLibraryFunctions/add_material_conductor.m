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

function add_material_conductor(MatName,library)

if nargin==1
    library = 'custom';
end

switch library
    case 'built-in'
        libName = 'materialLibrary\conductor_material.mat';
    case 'custom'
        libName = 'materialLibrary\custom_conductor.mat';
end


prompt={'Conductivity [S/m]',...
        'Density [kg*m^-3]',...
        'Thermal coefficient [1/°C]'};
name='New conductor material';
numlines=1;
defaultanswer={'580000000','8900','0.004'};

answer=inputdlg(prompt,name,numlines,defaultanswer);

mat.MatName   = MatName;
mat.sigma     = eval(answer{1});
mat.kgm3      = eval(answer{2});
mat.alpha     = eval(answer{3});

load(libName)
ind=length(MatList);
ind=ind+1;
MatList{ind}=MatName;
MatLib{ind}=mat;

button = questdlg('Save?','SELECT','Yes','No','Yes');
if isequal(button,'Yes')
    save(libName,'MatList','MatLib');
    disp(['Material added to ' library ' Material Library'])
else
    disp('material not saved')
end