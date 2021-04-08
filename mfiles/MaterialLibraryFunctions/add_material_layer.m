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

function add_material_layer(MatName,library)

if nargin==1
    library = 'custom';
end

switch library
    case 'built-in'
        libName = 'materialLibrary\layer_material.mat';
    case 'custom'
        libName = 'materialLibrary\custom_layer.mat';
end


prompt={'Relative permeability',...
        'Density [kg*m^-3]'...
        'Remanence [T]',...
        'Conductivity [S/m]',...
        'Reference temperature [°C]',...
        'Br=f(temperature) [T]',...
        'Bd=f(temperature) [T]'};    
name='New layer material';
numlines=1;
defaultanswer={'1','0','0','0','[20]','[0]','[0]'};
% defaultanswer={'1','0','0','0','0','[20]','[0]','[0]'};

answer=inputdlg(prompt,name,numlines,defaultanswer);

mat.MatName   = MatName;
mat.mu        = eval(answer{1});
mat.kgm3      = eval(answer{2});
mat.Br        = eval(answer{3});
%mat.Hc        = eval(answer{3});
mat.sigmaPM   = eval(answer{4});
mat.temp.temp = eval(answer{5});
mat.temp.Br   = eval(answer{6});
mat.temp.Bd   = eval(answer{7});
% mat.Hc        = eval(answer{4});
% mat.sigmaPM   = eval(answer{5});
% mat.temp.temp = eval(answer{6});
% mat.temp.Br   = eval(answer{7});
% mat.temp.Bd   = eval(answer{8});

mat.Hc = mat.Br/(mat.mu*(4e-7*pi));

load(libName)
ind=length(MatList);
ind=ind+1;
MatList{ind}=MatName;
MatLib{ind}=mat;

figure()
figSetting()
xlabel('$\theta$ [$^\circ$C]')
ylabel('[T]')
plot(mat.temp.temp,mat.temp.Br,'-go','DisplayName','$B_r$');
plot(mat.temp.temp,mat.temp.Bd,'-ro','DisplayName','$B_d$');
legend('show')

button = questdlg('Save?','SELECT','Yes','No','Yes');
if isequal(button,'Yes')
    save(libName,'MatList','MatLib');
    disp(['Material added to ' library ' Material Library'])
else
    disp('material not saved')
end