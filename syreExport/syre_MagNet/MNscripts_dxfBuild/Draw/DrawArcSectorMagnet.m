function [xC,thetaMean] = DrawArcSectorMagnet(magData,x0,rhoMin,rhoMax,thetaS,thetaE)

% DrawArcSectorMagnet.m [v1.00.00 (30-11-2012)]
% Draws an arc sector
% ====================================================================================
% Syntax: [xC,thetaMean] = DrawArcSectorMagnet(magData,x0,rhoMin,rhoMax,thetaS,thetaE)
% Input:
%          - magData: Magnet's data structure
%          - x0:      center of the arc
%          - rhoMin:  min radius
%          - rhoMax:  max radius
%          - thetaS:  starting angle (degrees)
%          - thetaE:  ending angle (degrees)
%
% Output:
%          - xC: center of the arc sector
%          - thetaMean: mean angle of the arc sector
% ====================================================================================

% degrees to radians
thetaSdeg = thetaS;
thetaEdeg = thetaE;
thetaS = thetaS*pi/180;
thetaE = thetaE*pi/180;

% polar to cartesian
pt1(1) = rhoMin*cos(thetaS);
pt1(2) = rhoMin*sin(thetaS);
pt2(1) = rhoMax*cos(thetaS);
pt2(2) = rhoMax*sin(thetaS);

pt3(1) = rhoMax*cos(thetaE);
pt3(2) = rhoMax*sin(thetaE);
pt4(1) = rhoMin*cos(thetaE);
pt4(2) = rhoMin*sin(thetaE);

% compute euclidean distance between pt1 and pt4
euDist = sqrt(abs(pt1(1)-pt4(1))^2+abs(pt1(2)-pt4(2))^2);
DrawLineMagnet(magData,pt1,pt2)
DrawArcMagnet(magData,x0,rhoMax,thetaSdeg,thetaEdeg);

if euDist>1e-6
    % Draw the lower sector
    DrawLineMagnet(magData,pt3,pt4);
    DrawArcMagnet(magData,x0,rhoMin,thetaSdeg,thetaEdeg);
else
    % shrink to a single point
    DrawLineMagnet(magData,pt3,pt1);
end

thetaMean = mean([thetaS thetaE]);
rhoMean = mean([rhoMin rhoMax]);
xC = [rhoMean*cos(thetaMean) rhoMean*sin(thetaMean)];
thetaMean = thetaMean*180/pi;