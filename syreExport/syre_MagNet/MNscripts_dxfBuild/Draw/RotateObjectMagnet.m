function RotateObjectMagnet(magData,objectList,rotAngle,center,rotAxis)

% RotateObjectMagnet.m [v1.00.00 (30-11-2012)]
% Rotates a list of objects
% ==================================================================================
% Syntax: RotateObjectMagnet(h,objectList,rotAngle,center,rotAxis)
% Input:
%          - magData: Magnet's data structure
%          - objectList: cell array containng the names of the objects to be rotated
%          - rotAngle: rotation angle
%          - center: center of rotation
%          - rotAxis: axis of rotation
% ==================================================================================


mh = magData.magnetHandler;
nObject = length(objectList);
objectListString = '';
for k = 1 : nObject
    objectListString = [objectListString,',"',objectList{k},'"'];
end
objectListString(1)=[]; % cancello la prima virgola
paramString = [num2str(center(1)),',',num2str(center(2)),',',num2str(center(3)),...
               ',',num2str(rotAxis(1)),',',num2str(rotAxis(2)),',',num2str(rotAxis(3)),...
               ',',num2str(rotAngle),',',num2str(1)];
invoke(mh,'processCommand',['CALL getDocument().rotateComponent(Array(',objectListString,'),',...
                           paramString,')']);
