% Copyright 2019
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [rotor,BLKLABELSrot]=drawPoleFBS(geo,mat,th_FBS)

geo.delta_FBS=th_FBS;
geo.dalpha(1)=geo.dalpha(1)+geo.delta_FBS/2*180/pi;
geo.x0=geo.r/cos(pi/2/geo.p+th_FBS);

mat.LayerMag.Br = mat.LayerMag.Br.*ones(1,geo.nlay);

if strcmp(geo.RotType,'Circular')
    [geo,mat,temp]=nodes_rotor_Circ_dx(geo,mat);
    rotor = build_matrix_Circ_dx(temp,geo);
elseif strcmp(geo.RotType,'Seg')
    [geo,mat,temp]=nodes_rotor_Seg(geo,mat);
    rotor=build_matrix_Seg(temp,geo);
else
    str=['Flux Barrier Shift not supported for ' geo.RotType ' geometry'];
    error(str);
end
    
% specchio il mezzo polo
rotNeg=rotor;
rotNeg(:,[2 4 6 7])=-rotor(:,[2 4 6 7]);
rotor=[rotor;rotNeg];

mesh_res = dimMesh(geo,'singt');

% find the centers of all blocks
BarCenter = defineBlockCenters(temp,mesh_res,geo);
% Assign label names
%BarName = defineBlockNames(temp,geo,mat);

codBound_periodic = -10; % simulate full machine

% shaft boundary
[xShaftBound1,yShaftBound1] = rot_point(mean([0,geo.Ar]),0,-90/geo.p*pi/180);
[xShaftBound2,yShaftBound2] = rot_point(mean([0,geo.Ar]),0,(geo.ps-1/2)*180/geo.p*pi/180);
% rotor boundary
[xRotBound1,yRotBound1] = rot_point(mean([geo.Ar,geo.r]),0,-90/geo.p*pi/180);
[xRotBound2,yRotBound2] = rot_point(mean([geo.Ar,geo.r]),0,(geo.ps-1/2)*180/geo.p*pi/180);

% Rotate rotor in zero position
rotor = rotateMatrix(rotor,pi/2/geo.p);

rotor = checkPlotMatrix(rotor,1e-9);

%%% Block centers %%%
BLKLABELSrot.xy     =   BarCenter;
%BLKLABELSrot.BarName =   BarName';
% Rotate block labels selection points
[xtemp,ytemp]=rot_point(BLKLABELSrot.xy(:,1),BLKLABELSrot.xy(:,2),90/geo.p*pi/180);
BLKLABELSrot.xy=[xtemp,ytemp,BLKLABELSrot.xy(:,3:end)];
clear xtemp ytemp;

% Magnetization direction rotation
xtemp=cos(atan2(BLKLABELSrot.xy(:,7),BLKLABELSrot.xy(:,6))+(pi/2/geo.p-eps));
ytemp=sin(atan2(BLKLABELSrot.xy(:,7),BLKLABELSrot.xy(:,6))+(pi/2/geo.p-eps));
BLKLABELSrot.xy(:,6)=xtemp;
BLKLABELSrot.xy(:,7)=ytemp;
clear xtemp ytemp;

%%% Boundaries %%%
BLKLABELSrot.boundary = [xShaftBound1,yShaftBound1,codBound_periodic;
    xShaftBound2,yShaftBound2,codBound_periodic;
    xRotBound1,yRotBound1,codBound_periodic;
    xRotBound2,yRotBound2,codBound_periodic];
% Rotate boundary selection points
[xtemp,ytemp]=rot_point(BLKLABELSrot.boundary(:,1),BLKLABELSrot.boundary(:,2),90/geo.p*pi/180);
BLKLABELSrot.boundary=[xtemp,ytemp,BLKLABELSrot.boundary(:,3:end)];


