% Copyright 2019
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

function FEMM_initialize(geo,mat)

% initial settings of FEMM problem

% problem definition
newdocument(0);
mi_probdef(0,'millimeters','planar',1e-8,geo.l,15);

% add iron in use
mi_addmaterial(mat.Stator.MatName);

for ii=1:length(mat.Stator.BH(:,1))
    mi_addbhpoint(mat.Stator.MatName,mat.Stator.BH(ii,1),mat.Stator.BH(ii,2));
end

if ~strcmp(mat.Stator.MatName,mat.Rotor.MatName)
    mi_addmaterial(mat.Rotor.MatName);
    for ii=1:length(mat.Rotor.BH(:,1))
        mi_addbhpoint(mat.Rotor.MatName,mat.Rotor.BH(ii,1),mat.Rotor.BH(ii,2));
    end
end

if ~strcmp(mat.Stator.MatName,mat.Shaft.MatName)
    mi_addmaterial(mat.Shaft.MatName);
    for ii=1:length(mat.Shaft.BH(:,1))
        mi_addbhpoint(mat.Shaft.MatName,mat.Shaft.BH(ii,1),mat.Shaft.BH(ii,2));
    end
end

% add conductor
mi_addmaterial(mat.SlotCond.MatName);
mi_modifymaterial(mat.SlotCond.MatName,5,mat.SlotCond.sigma/1e6);

% add barrier 
if isfield(mat.LayerMag,'BH')
    mi_addmaterial(mat.LayerMag.MatName);

    for ii=1:length(mat.LayerMag.BH(:,1))
        mi_addbhpoint(mat.LayerMag.MatName,mat.LayerMag.BH(ii,1),mat.LayerMag.BH(ii,2));
    end
    
else
    mi_addmaterial(mat.LayerMag.MatName,mat.LayerMag.mu,mat.LayerMag.mu,mat.LayerMag.Hc(1),0,mat.LayerMag.sigmaPM/1e6);
end

% add air
mi_addmaterial('Air',1,1,0,0);