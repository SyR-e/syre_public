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

function [mat] = material_properties_iron(MatName)

load('materialLibrary\iron_material.mat')

ind=0;

for ii=1:length(MatList)
    if strcmp(MatList{ii},MatName)
        ind=ii;
    end
end

matListBase = MatList;

if ind~=0
    mat=MatLib{ind};
else
    load('materialLibrary\custom_iron.mat')
    ind=0;
    for ii=1:length(MatList)
        if strcmp(MatList{ii},MatName)
            ind=ii;
        end
    end
    
    matListCustom = MatList;
    
    if ind~=0
        mat=MatLib{ind};
    else
        mat.MatName = MatName;
        mat.MatList = [matListBase matListCustom];
    end
end
