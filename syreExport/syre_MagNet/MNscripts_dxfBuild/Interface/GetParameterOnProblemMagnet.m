function p = GetParameterOnProblemMagnet(magData,ID,parName)
% GetParameterOnProblemMagnet.m [v1.00.00 (03-04-2013)]
% Gets the value of a parameter considering its instance on a specific
% problem ID.
% =========================================================================
% Syntax: p = GetParameterOnProblemMagnet(magData,pID,parName)
% Input:
%          - magData: Magnet's initial data structure
%          - ID:      Magnet's problem ID
%          - parName: name of the parameter to be read
%
% Output:
%          - p: parameter value
% =========================================================================

mh = magData.magnetHandler;
invoke(mh,'processCommand',['CALL getDocument().getProblem(',num2str(ID),').getParameter("","',parName,'",expression)']);
invoke(mh, 'processCommand', 'Call setVariant(0, expression)');
ps = invoke(mh, 'getVariant', 0);

% Returned could be a cell array or a string
if (ischar(ps))||(isnumeric(ps))
    % output is a string or a number
    p = ps;
elseif isa(ps,'cell')
    % output is a cell array of numeric values
    n = length(ps);
    p = zeros(1,n);
    for k = 1 : n
        p(k) = ps{k};
    end
end
