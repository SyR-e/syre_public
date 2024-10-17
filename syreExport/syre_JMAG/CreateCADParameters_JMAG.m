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

function CreateCADParameters_JMAG(CADEquationTable,VariableName, VariableType, VariableExpression, VariableDescription)
% 'CAD Function to create a CAD Parameters
    CADEquationTable.EditStart();
    CADEquationTable.AddEquation(VariableName);
    CADEquationTable.GetEquation(VariableName).SetType(VariableType);
    CADEquationTable.GetEquation(VariableName).SetExpression(VariableExpression);
    CADEquationTable.GetEquation(VariableName).SetDescription(VariableDescription);
    CADEquationTable.EditEnd();
end
