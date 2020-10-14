function SetBdryTangentialMagnet(magData,bdryName,bdryFaces)

% SetBdryTangentialMagnet.m [v1.00.00 (30-11-2012)]
% Sets tangential (Dirichlet) boundary conditions
% =========================================================================
% SetBdryTangentialMagnet(magData,bdryName,bdryFaces)
%
% Input:
%          - magData:   Magnet's data structure
%          - bdryName:  name of the boundary condition
%          - bdryFaces: cell matrix containing N rows and 2 columns. Each
%                       row represents a surface interested by the boundary
%                       condition. 
%                       Example:
%                       bdryFaces{i}{1} = object name (string)
%                       bdryFaces{i}{2} = object face (number)
%
% =========================================================================

% Get magnet's handler
mh = magData.magnetHandler;

% Get document's handler
dh = magData.magnetHandler;

% Number of faces
nFaces = length(bdryFaces);

% Selection of involved faces
invoke(mh, 'processCommand',['ReDim ArrayOfValues(',num2str(nFaces-1),')']);
for k = 1 : nFaces
    objName    = bdryFaces{k}{1};
    faceNumber = bdryFaces{k}{2};
%     SelectObjectFaceMagnet(h,objName,faceNumber,1);
    faceString = [objName,',Face#',num2str(faceNumber)];
    invoke(mh, 'processCommand',['ArrayOfValues(',num2str(k-1),')="',faceString,'"']);
end

% Set boundary condition
invoke(mh, 'processCommand',['CALL getDocument().createBoundaryCondition(ArrayOfValues,"',bdryName,'")']);
invoke(mh, 'processCommand',['CALL getDocument().setMagneticFluxTangential("',bdryName,'")']);

    
 
