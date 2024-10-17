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

function FEMCoil_Creation_JMAG (circuit, phaseLabel, position)
        circuit.CreateComponent('Coil', strcat(phaseLabel, '-Phase Coil'));
        circuit.CreateInstance(strcat(phaseLabel, '-Phase Coil'), position(1), position(2));
        circuit.GetComponent(strcat(phaseLabel, '-Phase Coil')).SetValue('Turn', 'coilTurn');
        circuit.GetComponent(strcat(phaseLabel, '-Phase Coil')).SetValue('Resistance', 'coilRes');
        circuit.GetComponent(strcat(phaseLabel, '-Phase Coil')).SetValue('LeakageInductance', 'coilLind');
end

