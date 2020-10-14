function SetTransientOptions (magData,tStart,tStep,tStop,sources,storageTime)

% SetTransientOptions.m [v1.00.00 (17-05-2013)]
% Sets the Transient option in MagNet
% ===============================================================================================
% Syntax: SetTransientOptions (magData,tStart,tStep,tStop,sources)
%
% Input:
%          - magData:   Magnet's data structure
%          - tStart:    Initial time of the simulation
%          - tStep:     Simulation step
%          - tStop:     Final time of the simulation
%          - sources:   Sources on at the beginning of the simulation. It
%                       must be 'Yes' or 'No'
% ===============================================================================================

% Get magnet's handler
mh = magData.magnetHandler;

invoke(mh,'processCommand','CALL getDocument().beginUndoGroup("Set Properties", true)');
invoke(mh,'processCommand',['CALL getDocument().setFixedIntervalTimeSteps(',num2str(tStart),', ',num2str(tStop),', ',num2str(tStep),')']);
invoke(mh,'processCommand','CALL getDocument().deleteTimeStepMaximumDelta()');

invoke(mh,'processCommand',['Call getDocument().setTimeStepStorageStartTime(',num2str(storageTime),')']); % Inizia a calcolare la soluzione da t = storageTime

invoke(mh,'processCommand',['CALL getDocument().setParameter("", "SourcesOnAtTransientStart", "',sources,'", infoStringParameter)']);
invoke(mh,'processCommand','CALL getDocument().endUndoGroup()');