function DrawArcMagnetperPunti(magData,x0,y0,x1,y1,x2,y2)

% DrawArcMagnet.m [v1.00.00 (30-11-2012)]
% Draws an arc
% =========================================================================
% Syntax: DrawArcMagnet(magData,x0,rho,thetaS,thetaE)
% Input:
%          - magData: Magnet's data structure
%          - x0:      center of the arc
%          - rho:     radius of the arc
%          - thetaS:  starting angle (degrees)
%          - thetaE:  ending angle (degrees)
% =========================================================================



vh = magData.viewHandler;
invoke(vh,'newArc',x0,y0,x1,y1,x2,y2);
