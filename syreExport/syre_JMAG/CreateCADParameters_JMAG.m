% 'CAD Function to create a CAD Parameters

function CreateCADParameters(CADEquationTable,VariableName, VariableType, VariableExpression, VariableDescription)
    CADEquationTable.EditStart();
    CADEquationTable.AddEquation(VariableName);
    CADEquationTable.GetEquation(VariableName).SetType(VariableType);
    CADEquationTable.GetEquation(VariableName).SetExpression(VariableExpression);
    CADEquationTable.GetEquation(VariableName).SetDescription(VariableDescription);
    CADEquationTable.EditEnd();
end
