function magData = OpenMagnet(visibilityFlag)

% OpenMagnet.m [v1.00.00 (30-11-2012)]
% Opens an ActiveX channel between Matlab and Magnet
% =========================================================================
% Syntax: magData = OpenMagnet(visibilityFlag)
% Input:
%          - visibilityFlag: if 1 Magnet is visible, if 0 it isn't.
% Output:
%          - magData: Magnet's initial data structure
%
% =========================================================================

% Open ActiveX channel
h = actxserver('Magnet.application');
set(h, 'Visible', visibilityFlag);

magData.magnetHandler = h;