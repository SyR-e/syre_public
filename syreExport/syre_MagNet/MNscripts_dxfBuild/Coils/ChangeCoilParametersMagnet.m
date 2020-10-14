function ChangeCoilParametersMagnet (magData,coilName,coilCurrent,frequency,phaseAngle)

% ChangeCoilParametersMagnet.m [v1.00.00 (20-05-2013)]
% This function change the phase angle parameter in the current expression
% =========================================================================================
% Syntax: ChangeCoilParametersMagnet (magData,coilName,coilCurrent,frequency,phaseAngle)
% Input:
%          - magData:     Magnet's data structure
%          - coilName:    name of the coil
%          - coilCurrent: coil current
%          - frequency:   electrical frequency of the waveform [Hz]
%          - phaseAngle:  current phase [electrical deg]
% =========================================================================

dh = magData.documentHandler;
mh = magData.magnetHandler;

invoke(mh,'processCommand',['CALL getDocument().getView().selectObject("',coilName,'", infoSetSelection)']);
invoke(mh,'processCommand',['CALL getDocument().beginUndoGroup("Set ',coilName,' Properties", true)']);
invoke(mh,'processCommand',['CALL getDocument().setParameter("',coilName,'", "WaveFormValues", "[0, ',num2str(coilCurrent),', ',num2str(frequency),', 0, 0, ',num2str(phaseAngle),']", infoArrayParameter)']);
invoke(mh,'processCommand',('CALL getDocument().endUndoGroup()'));


