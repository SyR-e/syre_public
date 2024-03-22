function createDesignParameters_JMAG(Design_parameters,ParameterName,type,expression,ParameterDescription)
    Design_parameters.AddEquation(ParameterName);
    Design_parameters.GetEquation(ParameterName).SetType(type);%type = 0, expression-dependent, type = 1, independent
if type == 1
    Design_parameters.GetEquation(ParameterName).SetExpression(expression);
else
    Design_parameters.GetEquation(ParameterName).SetExpression(strcat(num2str(expression)));
end
    Design_parameters.GetEquation(ParameterName).SetDescription (ParameterDescription);
end