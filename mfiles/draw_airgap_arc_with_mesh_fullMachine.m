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

% rebuilds the three arcs on the sliding surface between the stator and the
% rotor
% re-assigns the AP boundary conditions
% re-assigns the angular resolution for having a regular airgap mesh

function draw_airgap_arc_with_mesh_fullMachine(geo,th_m,fem)

mi_clearselected

p = geo.p;
r = geo.r;
g = geo.g;
q = geo.q;
% ns = geo.ns;
ps=geo.ps;
% Mezzo passo cava
pc = 360/(6*q*p)/2;

res=fem.res_traf;

% si fa riferimento ai gradi meccanici della porzione di rot simulata
gradi_da_sim = 180;

group = 20;

x0 = r+1/2*g;
y0 = 0;

res_deg=res;
angoli_bordo_mobile = [-pc,(-pc + gradi_da_sim)];

[x1,y1] = rot_point(x0,y0,angoli_bordo_mobile(1)*pi/180);
[x2,y2] = rot_point(x0,y0,angoli_bordo_mobile(2)*pi/180);

% num_seg=delta_angoli_bordo_mobile(2)*pi/180*x0*res/ps;
mi_addnode(x1,y1);mi_addnode(x2,y2);
mi_addarc(x1,y1,x2,y2,180,res);
mi_addarc(x2,y2,x1,y1,180,res);
[x,y] = rot_point(x2,y2,0.5*pi);
mi_selectarcsegment(x,y);
[x,y] = rot_point(x2,y2,-0.5*pi);
mi_selectarcsegment(x,y);
mi_setarcsegmentprop(res_deg, 'None', 0, group);
mi_clearselected

mi_addblocklabel(r+1/4*g,0);
mi_selectlabel(r+1/4*g,0);
mi_setblockprop('Air', 0, fem.res_traf, 'None', 0, 1, 1);
mi_clearselected;
mi_addblocklabel(r+3/4*g,0);
mi_selectlabel(r+3/4*g,0);
mi_setblockprop('Air', 0, fem.res_traf, 'None', 0, 1, 1);
mi_clearselected;

end



