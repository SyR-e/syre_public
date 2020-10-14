function magDataUpt = OpenDocumentMagnet(magData,docName)

% OpenDocumentMagnet.m [v1.01.00 (03-05-2013)]
% Opens a Magnet's document
% =========================================================================
% Syntax: magDataUpt = OpenDocumentMagnet(magData,docName)
% Input:
%         - magData:    Magnet's data structure
%         - docName:    path and name of the document to be open
% Output:
%         - magDataUpt: updated Magnet's data structure
% =========================================================================

h = magData.magnetHandler;
magDataUpt = magData;

% Open the document and get document handler
dh = invoke(h, 'openDocument', docName);

% Get current view
vh = invoke(dh, 'getCurrentView');

% Get Magnet's constants
Consts = invoke(h, 'getConstants');

% Update Magnet Data Strucure
magDataUpt.FileName = docName;
magDataUpt.documentHandler = dh;
magDataUpt.viewHandler = vh;
magDataUpt.magnetConstants = Consts;
magDataUpt.componentList = GetComponentsNameMagnet(magDataUpt);