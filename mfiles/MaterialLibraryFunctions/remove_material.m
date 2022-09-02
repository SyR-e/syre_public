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

function remove_material(MatName,type,library)
% 
% remove_material(MatName,type)
% 

if nargin==2
    library = 'custom';
end

switch library
    case 'built-in'
        switch type
            case 'Iron'
                filename='materialLibrary\iron_material.mat';
            case 'Conductor'
                 filename='materialLibrary\conductor_material.mat';
            case 'Layer'
                filename='materialLibrary\layer_material.mat';
            case 'Sleeve'
                filename='materialLibrary\sleeve_material.mat';
        end
    case 'custom'
       switch type
            case 'Iron'
                filename='materialLibrary\custom_iron.mat';
            case 'Conductor'
                 filename='materialLibrary\custom_conductor.mat';
            case 'Layer'
                filename='materialLibrary\custom_layer.mat';
           case 'Sleeve'
               filename='materialLibrary\custom_sleeve.mat';
        end
end

% if strcmp(type,'Iron')
%     filename='materialLibrary\iron_material.mat';
% elseif strcmp(type,'Conductor')
%     filename='materialLibrary\conductor_material.mat';
% elseif strcmp(type,'Layer')
%     filename='materialLibrary\layer_material.mat';
% else
%     disp('Wrong material type!')
%     disp('Allowed material type:')
%     disp('Iron')
%     disp('Conductor')
%     disp('Layer')
%     error('Insert a correct material type!')
% end

load(filename)

index=0;
for ii=1:length(MatList)
    if strcmp(MatName,MatList{ii})
        index=ii;
    end
end

if index==0
    disp([MatName ' not found in ' library ' ' type ' material library'])
    disp('Materials present:')
    for ii=1:length(MatList)
        disp(['- ' MatList{ii}]);
    end
    error('Insert a correct material name')
end

MatListOld=MatList;
MatLibOld=MatLib;
clear MatList MatLib

ind=1;

for ii=1:length(MatListOld)
    if ii~=index
        MatList{ind}=MatListOld{ii};
        MatLib{ind}=MatLibOld{ii};
        ind=ind+1;
    end
end

save(filename,'MatList','MatLib');

disp([MatName ' deleted from ' library ' ' type ' material library'])






