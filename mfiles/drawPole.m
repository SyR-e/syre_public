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

function [geo,temp,mat] = drawPole(geo,mat,fem)
% 
% This function draw a single standard pole, in zero position (as drawn in
% SyR-e), with labels matrix (BLKLABELSrot), 
% 

geo.delta_FBS=0; % no pole deformation
flagVtype = 1; % if 0, use Marco Gallo Vtype version, else use Simone Ferrari version (ready for syrmDesign)
if ~strcmp(geo.RotType,'SPM') || ~strcmp(geo.RotType,'Spoke-type')
    mat.LayerMag.Br = mat.LayerMag.Br.*ones(1,geo.nlay);   % replicate Br in case it is scalar
end

% 1) Find the design points (nodes_rotor_xxx), build rotor matrix
switch geo.RotType
    case 'Circular'
        % build nodes, lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_Circ_dx(geo,mat);
        rotor=build_matrix_Circ_dx(temp,geo);
    case 'ISeg'
        % build nodes, lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_ISeg(geo,mat);
        rotor=build_matrix_Seg(temp,geo);
    case 'Seg'
        % build nodes, lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_Seg(geo,mat);
        rotor=build_matrix_Seg(temp,geo);
    case 'Fluid'
        % build nodes, lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_Fluid(geo,mat);
        rotor=build_matrix_Fluid(temp,geo);
    case 'SPM'
        % build nodes, lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_SPM(geo,mat);
        rotor=build_matrix_SPM(temp,geo);
    case 'Vtype'
        % build nodes, lines and arcs for half a pole
        if flagVtype
            [geo,mat,temp] = nodes_rotor_Vtype_v3(geo,mat);
            rotor=build_matrix_Vtype_v3(temp,geo);
            %warning('Vtype geometry v2 under development!!!')
        else
            [geo,mat,temp] = nodes_rotor_Vtype(geo,mat);
            rotor=build_matrix_Vtype(temp,geo);
        end
    case 'Spoke-type'
        % build nodes , lines and arcs for half a pole
        [geo,mat,temp] = nodes_rotor_Spoke(geo,mat);
        rotor = build_matrix_Spoke(temp,geo);
end

% check about the PM area (to add to the nodes_rotor_xxx.m functions)
if ~isfield(geo,'AreaCMax')
    geo.AreaCMax=zeros(1,geo.nlay);
    geo.AreaEMax=zeros(1,geo.nlay);
end

if size(rotor,2)==7
    rotor=[rotor, ones(size(rotor,1),2)];
end

%2) mirror the half pole
rotNeg=rotor;
rotNeg(:,[2 4 6 7]) = -rotor(:,[2 4 6 7]);
rotNeg(:,9) = rotNeg(:,9)+max(rotor(:,9));
rotor=[rotor;rotNeg];

%3) find the centers of all blocks (PMs and air)
BarCenter = defineBlockCenters(temp,fem,geo);

% Rotate rotor in zero position (rotor matrix and labels)

rotor = rotateMatrix(rotor,90/geo.p*pi/180);
rotor = checkPlotMatrix(rotor,1e-9);

% Rotate block labels selection points
[xtemp,ytemp] = rot_point(BarCenter(:,1),BarCenter(:,2),90/geo.p*pi/180);
BarCenter     = [xtemp,ytemp,BarCenter(:,3:end)];
clear xtemp ytemp;

% Magnetization direction rotation
xtemp=cos(atan2(BarCenter(:,7),BarCenter(:,6))+(pi/2/geo.p-eps));
ytemp=sin(atan2(BarCenter(:,7),BarCenter(:,6))+(pi/2/geo.p-eps));
BarCenter(:,6)=xtemp;
BarCenter(:,7)=ytemp;
clear xtemp ytemp;

%%% OUTPUT DATA %%%
%%%%%%%%%%%%%%%%%%%

%%% Block centers %%%
geo.BLKLABELS.rotore.xy = BarCenter;
geo.rotor = rotor;











