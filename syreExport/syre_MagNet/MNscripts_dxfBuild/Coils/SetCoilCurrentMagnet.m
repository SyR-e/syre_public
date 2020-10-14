function SetCoilCurrentMagnet(magData,coilName,coilCurrent)

% SetCoilCurrentMagnet.m [v1.00.00 (30-11-2012)]
% Creates a simple coil
% =========================================================================================
% Syntax: SetCoilCurrentMagnet(magData,coilName,coilCurrent)
% Input:
%          - magData:     Magnet's data structure
%          - coilName:    name of the coil
%          - coilCurrent: current amplitude (magnetostatic)
% =========================================================================



dh = magData.documentHandler;

% Set coil properties
invoke(dh,'setCoilCurrent',coilName,coilCurrent,0);