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


function Coil_DirectionDefinition_JMAG (geo,study,phaseLabel,PhP,PhN)
% coil direction Definition: Upward or Downward
rowP = find(geo.win.avv == PhP);rowN = find(geo.win.avv == PhN);
% 'Create FEM Coil conditions for U-phase
phase_name = strcat(phaseLabel, '-Phase Coil');
 study.CreateCondition('FEMCoil', phase_name);
 FEM_condition = study.GetCondition(phase_name);
 FEM_condition.SetLink(phase_name);
 FEM_condition.SetName(phase_name);
for slot_ID = 1:length(rowP)
                 if ~isempty(rowP)
 %coil set definition
 FEM_condition.CreateSubCondition('FEMCoilData', strcat(phaseLabel,'p',num2str(slot_ID)));
 FEM_condition.GetSubCondition(strcat(phaseLabel,'p',num2str(slot_ID))).SetValue('Direction2D', 0);
 sel = FEM_condition.GetSubCondition(strcat(phaseLabel,'p',num2str(slot_ID))).GetSelection();
 sel.SelectPart(strcat(phaseLabel,'p',num2str(slot_ID)));
 FEM_condition.GetSubCondition(strcat(phaseLabel,'p',num2str(slot_ID))).AddSelected(sel);

                 end
end
for slot_ID = 1:length(rowN)
                 if ~isempty(rowN)
 FEM_condition.CreateSubCondition('FEMCoilData', strcat(phaseLabel,'n',num2str(slot_ID)));
 FEM_condition.GetSubCondition(strcat(phaseLabel,'n',num2str(slot_ID))).SetValue('Direction2D', 1);
 sel = FEM_condition.GetSubCondition(strcat(phaseLabel,'n',num2str(slot_ID))).GetSelection();
 sel.SelectPart(strcat(phaseLabel,'n',num2str(slot_ID)));
 FEM_condition.GetSubCondition(strcat(phaseLabel,'n',num2str(slot_ID))).AddSelected(sel);

                 end
end
FEM_condition.RemoveSubCondition(0);
end

