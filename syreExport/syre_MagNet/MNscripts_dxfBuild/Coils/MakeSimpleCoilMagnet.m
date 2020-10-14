function MakeSimpleCoilMagnet(magData,slotList,coilName)

% MakeSimpleCoilMagnet.m [v1.00.00 (30-11-2012)]
% Creates a simple coil
% =========================================================================================
% Syntax: MakeSimpleCoilMagnet(magData,slotList,coilName)
% Input:
%          - magData:   Magnet's data structure
%          - slotList:  cell array (2 component) containing the names of the objects which
%                       make the coil. The first object is the positive
%                       terminal while the second is the negative one.
%          - coilName:  name of the coil
% =========================================================================================

mh = magData.magnetHandler;

% Set coil's slots
invoke(mh,'processCommand','REDIM ArrayOfValues(1)');
invoke(mh,'processCommand',['ArrayOfValues(0)= "',slotList{1},'"']);
invoke(mh,'processCommand',['ArrayOfValues(1)= "',slotList{2},'"']);

% Create coil
invoke(mh,'processCommand','CALL getDocument().makeSimpleCoil(1, ArrayOfValues)');
dh = magData.documentHandler;
coilPath = invoke(dh,'getAllCoilPaths'); % recupera tutti i nomi delle coils
invoke(dh,'renameObject',coilPath{end},coilName); %coilPath{end} è il nome dell'ultima coil creata 

% Set coil properties
invoke(dh,'setCoilCurrent',coilName,0,0); % coil current magnitude and angle set to zero

