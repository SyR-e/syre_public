function compList = GetComponentsNameMagnet(h)

% GetComponentsNameMagnet.m [v1.01.00 (03-05-2013)]
% Gets the list of components
% =========================================================================
% Syntax: compList = GetComponentsNameMagnet(h)
% Input:
%         - h:          Magnet's data structure
% Output:
%         - compList:   cell array containing the name of each component
% =========================================================================


mh = h.magnetHandler;
invoke(mh,'processCommand','ReDim path(0)');
invoke(h.magnetHandler,'processCommand','path(0) = getDocument().getAllComponentPaths()');
invoke(mh,'processCommand','CALL setVariant(0,path)');
Pj = invoke(mh,'getVariant',0);

% Conversion from cell matrix to cell array
numOfComp = length(Pj{1});
compList = {};
for k = 1 : numOfComp
    compList{k} = Pj{1}{k};
end