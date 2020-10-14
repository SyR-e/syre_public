function magDataUpt = MakeComponentMagnet(magData,x0,componentName,depth,materialName,magType,magDir,meshSize)

% MakeComponentMagnet.m [v1.00.00 (30-11-2012)]
% Creates a component
% =============================================================================================================
% Syntax: magDataUpt = MakeComponentMagnet(magData,x0,componentName,depth,materialName,magType,magDir,meshSize)
% Input:
%          - magData: Magnet's data structure
%          - x0: position of the object from wich the component is made
%          - componentName: name of the component
%          - depth: dimension along Z-axis
%          - materialName: name of the material to be used
%          - magType: type of magnetization ('None' if the component isn't a permanent magnet)
%          - magDir: magnetization direction
%          - meshSize: if >0 it specifies the maximum mesh element size
% Output:
%          - magDataUpt: updated Magnet's data structure
% =============================================================================================================


% Retrive handlers
vh = magData.viewHandler;
mh = magData.magnetHandler;
ch = magData.magnetConstants;

% Select object surface
xc = x0(1);
yc = x0(2);
invoke(vh, 'setScaledToFit',1);
invoke(vh, 'selectAt', xc, yc, get(ch,'infoSetSelection'), get(ch,'infoSliceSurface'));

if ~isequal(magType,'None')
    % Change material name adding magnetization properties
    nameString = materialName;
    typeString = ['Type=',magType];
    directionString =  ['Direction=[',num2str(magDir(1)),',',num2str(magDir(2)),',',num2str(magDir(3)),']'];
    
    materialName = [nameString,';',typeString,';',directionString];
end
    
% Make component in a line without magnetization
invoke(mh, 'processCommand', 'ReDim ComponentName(0)');
invoke(mh, 'processCommand', ['ComponentName(0) = "',componentName,'"']);
invoke(mh, 'processCommand', ['CALL getDocument().getView().makeComponentInALine(',num2str(depth),...
    ', ComponentName, "Name=',materialName,'", infoMakeComponentUnionSurfaces Or infoMakeComponentRemoveVertices)']);

magDataUpt = magData;
nItems = length(magDataUpt.componentList);
magDataUpt.componentList{nItems+1} = componentName;

% Set mesh size
if meshSize>0
    SetMeshElementSizeMagnet(magData,componentName,meshSize);
end