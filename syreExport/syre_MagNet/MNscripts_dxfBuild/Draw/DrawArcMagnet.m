function DrawArcMagnet(magData,x0,rho,thetaS,thetaE)

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


% degrees to radians
thetaS = thetaS*pi/180;
thetaE = thetaE*pi/180;

% polar to cartesian
xS(1) = rho*cos(thetaS);
xS(2) = rho*sin(thetaS);
xE(1) = rho*cos(thetaE);
xE(2) = rho*sin(thetaE);

vh = magData.viewHandler;
invoke(vh,'newArc',x0(1),x0(2),xS(1),xS(2),xE(1),xE(2));
