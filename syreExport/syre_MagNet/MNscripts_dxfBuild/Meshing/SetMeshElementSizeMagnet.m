function SetMeshElementSizeMagnet(magData,objectName,maxSize)

% SetMeshElementSizeMagnet.m [v1.00.00 (30-11-2012)]
% Set the maximum mesh element size of a component
% =========================================================================
% Syntax: SetMeshElementSizeMagnet(magData,objectName,maxSize)
% Input:
%         - magData:    Magnet's data structure
%         - objectName: name of the object
%         - maxSize:    maximum mesh element size
% =========================================================================

dh = magData.documentHandler;
invoke(dh,'setMaxElementSize',objectName,maxSize);
