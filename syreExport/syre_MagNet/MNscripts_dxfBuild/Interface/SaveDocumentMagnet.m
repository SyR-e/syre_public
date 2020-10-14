function [magDataUpt,f] = SaveDocumentMagnet(magData,docName)

% SaveDocumentMagnet.m [v1.00.00 (30-11-2012)]
% Saves a new Magnet's document
% =========================================================================
% Syntax: f = SaveDocumentMagnet(h,docName)
% Input:
%          - h: current Magnet's instance (ActiveX handler)
%          - docName: document's name
% Output:
%          - true in case of successful save, false otherwise
% =========================================================================
% Stefano Ettorre
% stefano.ettorre@aviogroup.com

h = magData.magnetHandler;
f = invoke(h, 'saveDocument', docName );

magDataUpt = magData;
magDataUpt.FileName = docName;
