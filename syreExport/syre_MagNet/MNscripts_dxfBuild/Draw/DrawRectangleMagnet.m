function xC = DrawRectangleMagnet(magData,x0,width,height)

% DrawRectangleMagnet.m [v1.00.00 (30-11-2012)]
% Draws a rectangle
% =========================================================================
% Syntax: xC = DrawRectangleMagnet(magData,x0,width,height)
% Input:
%          - magData: Magnet's data structure
%          - x0: left-lower corner
%          - width: magnet's width
%          - heigth: magnets' height
% Output:
%          - xC: center of the rectangle
% =========================================================================


pt1 = x0; 
pt2 = pt1 + [width 0];
pt3 = pt2 + [0 height];
pt4 = pt3 + [-width 0];

DrawLineMagnet(magData,pt1,pt2);
DrawLineMagnet(magData,pt2,pt3);
DrawLineMagnet(magData,pt3,pt4);
DrawLineMagnet(magData,pt4,pt1);

xC = zeros(1,2);
xC(1) = mean([pt1(1) pt2(1)]);
xC(2) = mean([pt1(2) pt4(2)]);