% Copyright 2024
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

function SetCoilMaterial_JMAG(study,coil_mat,eddycurrent_c,Ph)
        study.SetMaterialByName(strcat('Coil',Ph), coil_mat);
        study.GetMaterial(strcat('Coil',Ph)).SetValue('Laminated', 0);
        study.GetMaterial(strcat('Coil',Ph)).SetValue('EddyCurrentCalculation', eddycurrent_c);
        study.GetMaterial(strcat('Coil',Ph)).SetValue('UserConductivityType', 0);%'0=use material resistivity 1=Electric conductivity 2=Electric resistivity 3=Temperature dependant conductivity 4=Temperature dependant resistivity
end
% study.GetMaterial(strcat(coil_groupname,'_U')).SetValue('UserResistivityValue', coil_conductivity);
