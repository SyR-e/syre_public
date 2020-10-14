function SetBdryRoundPeriodicMagnet(magData,bdryName,bdryType,bdryFaces,rotAxis,center,rotAngle)

% SetBdryRoundPeriodicMagnet.m [v1.00.00 (30-11-2012)]
% Sets (a)periodic boundary conditions with rotational transformation
% ===============================================================================================
% Syntax: SetBdryRoundPeriodicMagnet(magData,bdryName,bdryType,bdryFaces,rotAxis,center,rotAngle)
%
% Input:
%          - magData:   Magnet's data structure
%          - bdryName:  name of the boundary condition
%          - bdryType:  type of boundary condition ('Odd' or 'Even')
%          - bdryFaces: cell matrix containing N rows and 2 columns. Each
%                       row represents a surface interested by the boundary
%                       condition. 
%                       Example:
%                       bdryFaces{i}{1} = object name (string)
%                       bdryFaces{i}{2} = object face (number)
%          - rotAxis:   axis of rotation for the boundary transformation
%          - center:    center of rotation for the boundary transformation
%          - rotAngle:  angle of rotation (degrees) for the boundary transformation
% ===============================================================================================

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

invoke(mh, 'processCommand','ReDim RotationAxis(2)');
invoke(mh, 'processCommand',['RotationAxis(0)=',num2str(rotAxis(1))]);
invoke(mh, 'processCommand',['RotationAxis(1)=',num2str(rotAxis(2))]);
invoke(mh, 'processCommand',['RotationAxis(2)=',num2str(rotAxis(3))]);
invoke(mh, 'processCommand','ReDim Center(2)');
invoke(mh, 'processCommand',['Center(0)=',num2str(center(1))]);
invoke(mh, 'processCommand',['Center(1)=',num2str(center(2))]);
invoke(mh, 'processCommand',['Center(2)=',num2str(center(3))]);

if isequal(bdryType,'Even')
    invoke(mh, 'processCommand',['Call getDocument().setEvenPeriodic("',bdryName,'", 1,',num2str(rotAngle),', RotationAxis, Null, Null, Center)']);
elseif isequal(bdryType,'Odd')
    invoke(mh, 'processCommand',['Call getDocument().setOddPeriodic("',bdryName,'", 1,',num2str(rotAngle),', RotationAxis, Null, Null, Center)']);
end