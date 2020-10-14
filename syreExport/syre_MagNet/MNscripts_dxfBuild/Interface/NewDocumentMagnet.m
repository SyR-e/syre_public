function magDataUpt = NewDocumentMagnet(magData)

% NewDocumentMagnet.m [v1.00.00 (30-11-2012)]
% Creates a new Magnet's document
% =========================================================================
% Syntax: dh = NewDocumentMagnet(h)
% Input:
%          - magData: current Magnet's data structure
% Output:
%          - magData: updated Magnet's data structure
% =========================================================================

h = magData.magnetHandler;

% New Document
dh = invoke(h, 'newDocument') ;

% Get current view
vh = invoke(dh, 'getCurrentView');

% Get Magnet's constants
Consts = invoke(h, 'getConstants');

% Set millimeters as default length unit
invoke(dh,'setDefaultLengthUnit','Millimeters');

magDataUpt = magData;
magDataUpt.documentHandler = dh;
magDataUpt.viewHandler = vh;
magDataUpt.magnetConstants = Consts;
magDataUpt.componentList = {};