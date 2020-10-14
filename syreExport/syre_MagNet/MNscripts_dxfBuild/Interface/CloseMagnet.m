function CloseMagnet(magData)

% CloseMagnet.m [v1.00.00 (30-11-2012)]
% Closes a Magnet's ActiveX channel
% =========================================================================
% Syntax: CloseMagnet(h)
% Input:
%          - magData: Magnet's data structure
%
% =========================================================================

% Close ActiveX channel
invoke(magData.magnetHandler, 'Exit');
release(magData.magnetHandler);
