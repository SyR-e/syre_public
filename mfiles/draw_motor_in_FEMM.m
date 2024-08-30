% Copyright 2014
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

function [geo,mat] = draw_motor_in_FEMM(geo,mat,pathname,filename)

% draw_motor_in_FEMM.m
% builds the .fem motor model with the rotor in position zero

% input: see function argument
% output:
% - updated geo (individual for the machine)
% - mot0.fem (th_m = 0, i123 = [0,0,0])

% % Sliding airgap definition
% geo.slidingGap=1;

openfemm(1);

FEMM_initialize(geo,mat);
% mi_probdef(0,'millimeters','planar',1e-8,geo.l,15);

eval_type = 'singt';    % default 
fem = dimMesh(geo,eval_type);

% calc winding factor (kavv) and rotor offset (phase1_offset)
[~,phase1_offset] = calcKwTh0(geo);
phase1_offset = phase1_offset+360/(6*geo.p*geo.q*geo.win.n3phase)/2*geo.p;    %first slot in 360/(6pq)/2 position

if strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype') || strcmp(geo.RotType,'SPM-Halbach')
    geo.axisType = 'PM';
    phase1_offset = phase1_offset-90;
elseif strcmp(geo.RotType,'Spoke-type')
    geo.axisType = 'PM';
    phase1_offset = phase1_offset;
else
    geo.axisType = 'SR';
    phase1_offset = phase1_offset;
end


% if strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype')
%     if strcmp(geo.axisType,'PM')
%         phase1_offset = phase1_offset - 90;   % valid for d axis on PM direction
%     else
%         phase1_offset = phase1_offset - 180;  % valid for -q axis on PM direction
%     end
% else
%     if strcmp(geo.axisType,'PM')
%         phase1_offset = phase1_offset + 90;   % valid for d axis on PM direction
%     else
%         phase1_offset = phase1_offset;
%     end
% end

%%% Axis type da A ad a - commenta gamma+90 - verifica spm to sr - 

geo.th0 = - phase1_offset;  % d- to alpha-axis offset [elt deg]

% if FBS, the mean q-axis must be aligned - used equation alignment with mean q-axis
alphaQs=[pi/(2*geo.p)];  % s=symmetric
alphaQa=[pi/(2*geo.p)];  % a=asymmetric
delta_FBS=geo.th_FBS*[-1 1];
for ii=2:length(delta_FBS)
    alphaQs(ii)=alphaQs(ii-1)+pi/geo.p;
    alphaQa(ii)=alphaQa(ii-1)+delta_FBS(ii-1)/2+delta_FBS(ii)/2+pi/geo.p;
end
alphaQas=alphaQa-alphaQs;
th_rot=geo.th_FBS/2;
th_rot=+th_rot-mean(alphaQas);

geo.th0=geo.th0-th_rot*geo.p*180/pi;

% boundary conditions: definition (assigned to segments later)
mi_addboundprop('A=0', 0, 0, 0, 0, 0, 0, 0, 0, 0);  % inner and outer circles
% Periodicity or Anti-Periodicity (2 x rotor + 2 x stator + 3 x airgap + 1 x APmove) APmove is the sliding contour
mi_addboundprop('APr1', 0, 0, 0, 0, 0, 0, 0, 0, geo.periodicity);   % shaft
mi_addboundprop('APr2', 0, 0, 0, 0, 0, 0, 0, 0, geo.periodicity);   % rotor iron
mi_addboundprop('APr3', 0, 0, 0, 0, 0, 0, 0, 0, geo.periodicity);   % sleeve
mi_addboundprop('APg1', 0, 0, 0, 0, 0, 0, 0, 0, geo.periodicity);
mi_addboundprop('APg2', 0, 0, 0, 0, 0, 0, 0, 0, geo.periodicity);
mi_addboundprop('APg3', 0, 0, 0, 0, 0, 0, 0, 0, geo.periodicity);
%mi_addboundprop('APmove', 0, 0, 0, 0, 0, 0, 0, 0, geo.periodicity);
mi_addboundprop('APs1', 0, 0, 0, 0, 0, 0, 0, 0, geo.periodicity);
mi_addboundprop('APs2', 0, 0, 0, 0, 0, 0, 0, 0, geo.periodicity);
mi_addboundprop('APs3', 0, 0, 0, 0, 0, 0, 0, 0, geo.periodicity);
mi_addboundprop('AGap', 0, 0, 0, 0, 0, 0, 0, 0, geo.periodicity+2,0,0);

% nodes
geo.x0 = geo.r/cos(pi/2/geo.p);
% geo.x0 = (geo.r-geo.hs)/cos(pi/2/geo.p);


[rotor,BLKLABELSrot,geo,mat] = ROTmatr(geo,fem,mat); % rotor and BLKLABELSrot describe the rotor
geo.rotor = rotor;

[geo,stator,BLKLABELSstat] = STATmatr(geo,fem); % statore and BLKLABELSstat describe the stator
geo.stator=stator;
% draw lines and arcs
draw_lines_arcs(rotor,2,fem.res);
draw_lines_arcs(stator,1,fem.res);

% check rotor arcs resolution (error for no ribs design and multiple arcs
% drawn)
% FEMM_checkArcResolution(rotor,2,fem.res)

% assign block labels
BLKLABELS.materials=geo.BLKLABELSmaterials;
BLKLABELS.rotore = BLKLABELSrot;
BLKLABELS.statore= BLKLABELSstat;
geo.BLKLABELS=BLKLABELS;
assign_block_prop_rot(BLKLABELS,geo,mat,fem,2);
assign_block_prop_stat(BLKLABELS,geo,fem,1);

% assign boundary conditions to the rotor
for ii=1:2
    mi_selectsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2));
    if (BLKLABELSrot.boundary(ii,3)==10)
        mi_setsegmentprop('APr1', 0, 1, 0, 2);
        mi_clearselected;
    elseif(BLKLABELSrot.boundary(ii,3)==0)
        mi_selectarcsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2))
        mi_setarcsegmentprop(fem.res, 'A=0', 0, 2);
        mi_clearselected;
    end
end

for ii=3:4
    mi_selectsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2));
    if (BLKLABELSrot.boundary(ii,3)==10)
        mi_setsegmentprop('APr2', 0, 1, 0, 2);
        mi_clearselected;
    elseif(BLKLABELSrot.boundary(ii,3)==0)
        mi_selectarcsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2))
        mi_setarcsegmentprop(fem.res, 'A=0', 0, 2);
        mi_clearselected;
    end
end

for ii=5:6
    mi_selectsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2));
    if (BLKLABELSrot.boundary(ii,3)==10)
        mi_setsegmentprop('APr3', 0, 1, 0, 2);
        mi_clearselected;
    elseif(BLKLABELSrot.boundary(ii,3)==0)
        mi_selectarcsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2))
        mi_setarcsegmentprop(fem.res, 'A=0', 0, 2);
        mi_clearselected;
    end
end

% assign boundary conditions to the stator
for ii=1:2
    mi_selectsegment(BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2));
    if (BLKLABELSstat.boundary(ii,3)==10)
        mi_setsegmentprop('APs1', 0, 1, 0, 1);
        mi_clearselected;
    elseif(BLKLABELSstat.boundary(ii,3)==0)
        mi_selectarcsegment(BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2));
        mi_setarcsegmentprop(fem.res, 'A=0', 0, 1);
        mi_clearselected;
    end
end

for ii=3:size(BLKLABELSstat.boundary,1)
    mi_selectsegment(BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2));
    if (BLKLABELSstat.boundary(ii,3)==10)
        mi_setsegmentprop('APs3', 0, 1, 0, 1);
        mi_clearselected;
    elseif(BLKLABELSstat.boundary(ii,3)==0)
        mi_selectarcsegment(BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2));
        mi_setarcsegmentprop(fem.res, 'A=0', 0, 1);
        mi_clearselected;
    end
end

% build the airgap (group 20)

draw_airgap(geo,fem);

% th_m0 = 0; % rotor aligned with the horizon
% if (geo.ps<2*geo.p)
%     AirGapBuild(geo.Qs,geo.ps,geo.p,geo.g,360/(6*geo.q*geo.p*geo.win.n3phase)/2,geo.r,fem.res_traf,1,2,geo.lm,BLKLABELSrot,geo.RotType);
% else
%     draw_airgap_arc_with_mesh_fullMachine(geo,th_m0,fem);
% end

geo.mesh_res=fem;
if isoctave() %OCT
    tmp = mat.LayerMag.Br;
    [~,ind,~] = unique(tmp,'first');
    mat.LayerMag.Br = tmp(sort(ind));
    %mat.LayerMag.Br = unique_oct(mat.LayerMag.Br,'stable');     
else    
    mat.LayerMag.Br = unique(mat.LayerMag.Br,'stable');   % reset geo.Br to input value (either scalar or size = nlay)
end

if ~isfield(geo,'custom')||(geo.custom==0)
%    filename = strrep(filename,'.fem','Syre.fem');
     mi_saveas([pathname filename]);
end
% mi_saveas([pathname strrep(filename,'.mat','.fem')]);
% mi_saveas([pathname filename]);