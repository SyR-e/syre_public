function DrawLineMagnet(magData,x1,x2)

% DrawLineMagnet.m [v1.00.00 (30-11-2012)]
% Draws a line
% =========================================================================
% Syntax: DrawLineMagnet(magData,x1,x2)
% Input:
%          - magData: Magnet's data structure
%          - x1:      starting point
%          - x2:      ending point
% =========================================================================

vh = magData.viewHandler;

invoke(vh, 'NewLine', x1(1),x1(2),x2(1),x2(2));
