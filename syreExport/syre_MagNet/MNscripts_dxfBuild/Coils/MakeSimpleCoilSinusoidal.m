function MakeSimpleCoilSinusoidal(magData,slotList,terminalList,coilName,type,coilCurrent,nTurns,frequency,phase)

% MakeSimpleCoilSinusoidal.m [v1.00.00 (17-05-2013)]
% Creates a simple coil
% =========================================================================================
% Syntax: MakeSimpleCoilSinusoidal(magData,slotList,terminalList,coilName,type,coilCurrent,nTurns,frequency,phase)

% Input:
%          - magData:       Magnet's data structure.
%          - slotList:      cell array (2 component) containing the names of the objects which
%                           make the coil.
%          - terminalList:  cell array (2 component) containing the names of the start and end  
%                           terminal faces of the coil.
%          - coilName:      name of the coil.
%          - type:          the coil type can be defined as "Voltage" or "Current" driven. The
%                           firsth letter of the word MUST BE CAPITAL.
%          - coilCurrent:   peak phase current in the coil.
%          - nTurns:        number of turns of the coil.
%          - frequency:     current electrical frequency in Hz.
%          - phase:         phase of the current in electrical degrees                 
%
% =========================================================================================




dh = magData.documentHandler;
mh = magData.magnetHandler;

% Set coil's slots
invoke(mh,'processCommand','REDIM ArrayOfValues(1)');
invoke(mh,'processCommand',['ArrayOfValues(0)= "',slotList{1},',',terminalList{1},'"']);
invoke(mh,'processCommand',['ArrayOfValues(1)= "',slotList{2},',',terminalList{2},'"']);

invoke(mh,'processCommand','CALL getDocument().makeSimpleCoil(1, ArrayOfValues)');
coilPath = invoke(dh,'getAllCoilPaths'); % recupera tutti i nomi delle coils
invoke(dh,'renameObject',coilPath{end},coilName); %coilPath{end} è il nome dell'ultima coil creata 


invoke(mh,'processCommand',['CALL getDocument().beginUndoGroup("Set ',coilName,' Properties", true)']);
invoke(mh,'processCommand',['CALL getDocument().setCoilSourceType("',coilName,'", info',type,'Driven)']);
invoke(mh,'processCommand',['CALL getDocument().setParameter("',coilName,'", "NumberOfTurns", "',num2str(nTurns),'", infoNumberParameter)']);
invoke(mh,'processCommand',['CALL getDocument().setParameter("',coilName,'", "WaveFormType", "SIN", infoStringParameter)']);
invoke(mh,'processCommand',['CALL getDocument().setParameter("',coilName,'", "WaveFormValues", "[0, ',num2str(coilCurrent),', ',num2str(frequency),', 0, 0, ',num2str(phase),']", infoArrayParameter)']);
invoke(mh,'processCommand','CALL getDocument().endUndoGroup()');
