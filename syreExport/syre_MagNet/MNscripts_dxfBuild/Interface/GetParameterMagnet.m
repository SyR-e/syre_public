function p = GetParameterMagnet(magData,parName)
% GetParameterMagnet.m [v1.00.00 (03-04-2013)]
% Gets the value of a parameter. The returned value IS NOT related to a
% specific problem ID. Should be used in preprocessing only.
% =========================================================================
% Syntax: p = GetParameterOnProblemMagnet(magData,pID,parName)
% Input:
%          - magData: Magnet's initial data structure
%          - parName: name of the parameter to be read
%
% Output:
%          - p: parameter value
% =========================================================================

mh = magData.magnetHandler;
invoke(mh,'processCommand',['CALL getDocument().getParameter("","',parName,'",expression)'])
invoke(mh, 'processCommand', 'Call setVariant(0, expression)');
ps = invoke(mh, 'getVariant', 0);

% Convert string to a number (if necessary)
p = str2num(ps);
if isempty(p)
    % returned value is a text and not a number
    p = ps;
end