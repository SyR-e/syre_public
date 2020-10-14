function SetParameterMagnet(magData,parName,parValue)
% SetParameterMagnet.m [v1.00.00 (03-04-2013)]
% Sets a parameter (numeric scalar or array and string)
% =========================================================================
% Syntax: SetParameterMagnet(magData,parName,parValue)
% Input:
%          - magData:  Magnet's initial data structure
%          - parName:  name of the parameter to be set
%          - parValue: value of the parameter to be set (could be a scalar,
%                      an array or a string)
%
% =========================================================================
mh = magData.magnetHandler;

% Check parameter type
chk = length(parValue)*isnumeric(parValue);

if chk==1
    % scalar numeric parameter (isnumeric = 1 and length = 1)
    pType = 'infoNumberParameter';
    pValue = num2str(parValue);
    
elseif chk>1
    % array numeric parameter (isnumeric = 1 and length > 1)
    pType = 'infoArrayParameter';
    arraySize = length(parValue);
    openingString = '[';
    closingString = ']';
    arrayString = [];
    for k = 1 : arraySize
        arrayString = [arrayString num2str(parValue(k)) ' '];
    end
    pValue = [openingString arrayString closingString];
    pValue(end-1) = [];
    
elseif chk==0
    % string parameter (isnumeric = 0)
    pType = 'infoStringParameter';
    pValue = parValue;
end

valueString = ['"',parName,'","',pValue,'",',pType];
cmdString = ['CALL getDocument().setParameter("",',valueString,')'];
invoke(mh,'processCommand',cmdString);