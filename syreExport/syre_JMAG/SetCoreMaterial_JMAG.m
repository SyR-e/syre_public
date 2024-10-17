%  Copyright 2024
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

function SetCoreMaterial_JMAG(study,core_groupname,core_mat,laminated_type,laminated_factor,hysteresisloop_index,eddycurrent_index)

    study.SetMaterialByName(core_groupname, core_mat);
    study.GetMaterial(core_groupname).SetValue('Laminated', laminated_type);
    study.GetMaterial(core_groupname).SetValue('LaminationFactor', laminated_factor);
    study.GetMaterial(core_groupname).SetValue('MagnetizationCorrection', 100);
if hysteresisloop_index==1
    study.GetMaterial(core_groupname).SetValue('UseMaterialHysteresisLoop', hysteresisloop_index);
end
    study.GetMaterial(core_groupname).SetValue('EddyCurrentCalculation', eddycurrent_index);
if eddycurrent_index==0
    study.GetMaterial(core_groupname).SetValue('ConductivityType', 0);
    study.GetMaterial(core_groupname).SetValue('UserConductivityType', 0);
end
end

