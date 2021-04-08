% Copyright 2019
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

function update_material_library(pathname,filename)
%
% update_material_library(pathname,filename)
% 

if nargin<2
    [filename,pathname,~]=uigetfile([cd '\.mat'],'Select a motor model');
end

flagMod=zeros(1,5);

% disp('Material update from:')
% disp([pathname filename])

mot=load([pathname filename]);

% check rotor iron
matTmp=material_properties_iron(mot.dataSet.RotorMaterial);
if ~isfield(matTmp,'kgm3')
    mat.MatName   = mot.mat.Rotor.MatName;
    mat.sigma_max = mot.mat.Rotor.sigma_max;
    mat.kgm3      = mot.mat.Rotor.kgm3;
    mat.alpha     = mot.mat.Rotor.alpha;
    mat.beta      = mot.mat.Rotor.beta;
    mat.kh        = mot.mat.Rotor.kh;
    mat.ke        = mot.mat.Rotor.ke;
    mat.BH        = mot.mat.Rotor.BH;
    if isfield(mot.mat.Rotor,'E')
        mat.E     = mot.mat.Rotor.E;
    else
        mat.E     = [];
    end
    
    %load('materialLibrary\iron_material.mat');
    load('materialLibrary\custom_iron.mat');
    ind=length(MatList);
    ind=ind+1;
    MatList{ind}=mat.MatName;
    MatLib{ind}=mat;
    
    save('materialLibrary\custom_iron.mat','MatList','MatLib');
    flagMod(2)=1;
    clear mat
end

% check stator iron
matTmp=material_properties_iron(mot.dataSet.StatorMaterial);
if ~isfield(matTmp,'kgm3')
    mat.MatName   = mot.mat.Stator.MatName;
    mat.sigma_max = mot.mat.Stator.sigma_max;
    mat.kgm3      = mot.mat.Stator.kgm3;
    mat.alpha     = mot.mat.Stator.alpha;
    mat.beta      = mot.mat.Stator.beta;
    mat.kh        = mot.mat.Stator.kh;
    mat.ke        = mot.mat.Stator.ke;
    mat.BH        = mot.mat.Stator.BH;
    if isfield(mot.mat.Stator,'E')
        mat.E     = mot.mat.Stator.E;
    else
        mat.E     = [];
    end
    
    load('materialLibrary\custom_iron.mat');
    ind=length(MatList);
    ind=ind+1;
    MatList{ind}=mat.MatName;
    MatLib{ind}=mat;
    
    save('materialLibrary\custom_iron.mat','MatList','MatLib');
    flagMod(1)=1;
    clear mat
end

% check shaft iron
matTmp=material_properties_iron(mot.dataSet.ShaftMaterial);
if (~isfield(matTmp,'kgm3')&&~strcmp(mot.dataSet.ShaftMaterial,'Air'))
    mat.MatName   = mot.mat.Shaft.MatName;
    mat.sigma_max = mot.mat.Shaft.sigma_max;
    mat.kgm3      = mot.mat.Shaft.kgm3;
    mat.alpha     = mot.mat.Shaft.alpha;
    mat.beta      = mot.mat.Shaft.beta;
    mat.kh        = mot.mat.Shaft.kh;
    mat.ke        = mot.mat.Shaft.ke;
    mat.BH        = mot.mat.Shaft.BH;
    
    load('materialLibrary\custom_iron.mat');
    ind=length(MatList);
    ind=ind+1;
    MatList{ind}=mat.MatName;
    MatLib{ind}=mat;
    
    save('materialLibrary\custom_iron.mat','MatList','MatLib');
    flagMod(3)=1;
    clear mat
end


% check PM material
matTmp=material_properties_layer(mot.dataSet.FluxBarrierMaterial);
if ~isfield(matTmp,'kgm3')
    mat.MatName = mot.mat.LayerMag.MatName;
    mat.mu      = mot.mat.LayerMag.mu;
    mat.kgm3    = mot.mat.LayerMag.kgm3;
    mat.Br      = mot.mat.LayerMag.Br;
    mat.sigmaPM = mot.mat.LayerMag.sigmaPM;
    if isfield(mot.mat.LayerMag,'temp')
        mat.temp    = mot.mat.LayerMag.temp;
    end
    mat.Hc      = mot.mat.LayerMag.Hc;
    
    load('materialLibrary\custom_layer.mat');
    ind=length(MatList);
    ind=ind+1;
    MatList{ind}=mat.MatName;
    MatLib{ind}=mat;
    
    save('materialLibrary\custom_layer.mat','MatList','MatLib');
    flagMod(4)=1;
    clear mat
end

% check slot material
matTmp=material_properties_conductor(mot.dataSet.SlotMaterial);
if ~isfield(matTmp,'kgm3')
    mat.MatName = mot.mat.SlotCond.MatName;
    mat.sigma   = mot.mat.SlotCond.sigma;
    mat.kgm3    = mot.mat.SlotCond.kgm3;
    mat.alpha   = mot.mat.SlotCond.alpha;
    
    load('materialLibrary\custom_conductor.mat');
    ind=length(MatList);
    ind=ind+1;
    MatList{ind}=mat.MatName;
    MatLib{ind}=mat;
    
    save('materialLibrary\custom_conductor.mat','MatList','MatLib');
    flagMod(5)=1;
    clear mat
end

if sum(flagMod)
    disp('Materials added to the custom material library:')
    for ii=1:length(flagMod)
        if flagMod(ii)
            switch ii
                case 1
                    disp(['- ' mot.dataSet.StatorMaterial '(iron)'])
                case 2
                    disp(['- ' mot.dataSet.RotorMaterial '(iron)'])
                case 3
                    disp(['- ' mot.dataSet.ShaftMaterial '(iron)'])
                case 4
                    disp(['- ' mot.dataSet.FluxBarrierMaterial '(PM)'])
                case 5
                    disp(['- ' mot.dataSet.SlotMaterial '(conductor)'])
            end
        end
    end
    beep
end







