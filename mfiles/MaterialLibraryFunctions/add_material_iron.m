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

function add_material_iron(MatName,library)

if nargin==1
    library = 'custom';
end

switch library
    case 'built-in'
        libName = 'materialLibrary\iron_material.mat';
    case 'custom'
        libName = 'materialLibrary\custom_iron.mat';
end

prompt={'Yield strength [MPa]',...
        'Density [kg*m^-3]',...
        'alpha (iron loss coeff)',...
        'beta (iron loss coeff)',...
        'kh (iron loss coeff)',...
        'ke (iron loss coeff)',...
        'Young Module [GPa]'};
name='New iron material';
numlines=1;
defaultanswer={'200','7800','0','0','0','0','200'};

answer=inputdlg(prompt,name,numlines,defaultanswer);

[filename,pathname,~]=uigetfile([cd '\.m'],'Load BH curve');
run([pathname filename]);
figure
plot(BH(:,2),BH(:,1));
if (BH(1,1)~=0 || BH(1,2)~=0)
    error('First point of BH curve must be (0,0)')
end

mat.MatName   = MatName;
mat.sigma_max = eval(answer{1});
mat.kgm3      = eval(answer{2});
mat.alpha     = eval(answer{3});
mat.beta      = eval(answer{4});
mat.kh        = eval(answer{5});
mat.ke        = eval(answer{6});
mat.BH        = BH;
mat.E         = eval(answer{7});

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