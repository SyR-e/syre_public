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

%% 2018/12/10: update with the new airgap boundaty condition (FEMM 28Feb2018)

function AirGapBuild(Qs,ps,p,g,pc,r,res_traf,groupStat,groupRot,lm,BLKLABELSrot,RotType)

groupStat=20;
groupRot=20;

% stator
xArcStat1 = r+2/3*g;
yArcStat1 = 0;
xArcStat2 = (r+2/3*g)*cos(pi/p*ps);
yArcStat2 = (r+2/3*g)*sin(pi/p*ps);

mi_drawarc(xArcStat1,yArcStat1,xArcStat2,yArcStat2,2*pc*Qs,res_traf);
mi_selectarcsegment(xArcStat1,yArcStat1);
mi_setarcsegmentprop(res_traf, 'AGap', 0, groupStat);
mi_selectnode(xArcStat1,yArcStat1); mi_setnodeprop('None',groupStat);
mi_selectnode(xArcStat2,yArcStat2); mi_setnodeprop('None',groupStat);
mi_clearselected

% linee al traferro arco traf --> statore
xTrafStat1 = r+g;
yTrafStat1 = 0;
xTrafStat2 = (r+g)*cos(pi/p*ps);
yTrafStat2 = (r+g)*sin(pi/p*ps);

mi_drawline(xArcStat1,yArcStat1,xTrafStat1,yTrafStat1);
mi_selectsegment(mean([xArcStat1 xTrafStat1]),mean([yArcStat1 yTrafStat1]));
mi_setsegmentprop('APg3', res_traf, 0, 0, groupStat);
mi_selectnode(xTrafStat1,yTrafStat1); mi_setnodeprop('None',groupStat);
mi_selectnode(xArcStat1,yArcStat1); mi_setnodeprop('None',groupStat);
mi_clearselected;

mi_drawline(xArcStat2,yArcStat2,xTrafStat2,yTrafStat2);
mi_selectsegment(mean([xArcStat2 xTrafStat2]),mean([yArcStat2 yTrafStat2]));
mi_setsegmentprop('APg3', res_traf, 0, 0, groupStat);
mi_selectnode(xTrafStat2,yTrafStat2); mi_setnodeprop('None',groupStat);
mi_selectnode(xArcStat2,yArcStat2); mi_setnodeprop('None',groupStat);
mi_clearselected;

[xAirTrafSt,yAirTrafSt] = rot_point(r+5/6*g,0,0.5*ps*180/p*pi/180);
mi_addblocklabel(xAirTrafSt,yAirTrafSt);
mi_selectlabel(xAirTrafSt,yAirTrafSt);
mi_setblockprop('Air', 0, res_traf, 'None', 0, groupStat, 1);
mi_clearselected;

% rotor
xArcRot1 = r+1/3*g;
yArcRot1 = 0;
xArcRot2 = (r+1/3*g)*cos(pi/p*ps);
yArcRot2 = (r+1/3*g)*sin(pi/p*ps);
if strcmp(RotType,'SPM')
    xTrafRot1 = r-lm;
    yTrafRot1 = 0;
    xTrafRot2 = (r-lm)*cos(pi/p*ps);
    yTrafRot2 = (r-lm)*sin(pi/p*ps);
else
    xTrafRot1 = r;
    yTrafRot1 = 0;
    xTrafRot2 = r*cos(pi/p*ps);
    yTrafRot2 = r*sin(pi/p*ps);
end

mi_drawarc(xArcRot1,yArcRot1,xArcRot2,yArcRot2,2*pc*Qs,res_traf);
mi_selectarcsegment(xArcRot1,yArcRot1);
mi_setarcsegmentprop(res_traf, 'AGap', 0, groupRot);
mi_selectnode(xArcRot1,yArcRot1); mi_setnodeprop('None',groupRot);
mi_selectnode(xArcRot2,yArcRot2); mi_setnodeprop('None',groupRot);
mi_clearselected

mi_drawline(xArcRot1,yArcRot1,xTrafRot1,yTrafRot1);
mi_selectsegment(mean([xArcRot1 xTrafRot1]),mean([yArcRot1 yTrafRot1]));
mi_setsegmentprop('APg1', res_traf, 0, 0, groupRot);
mi_selectnode(xTrafRot1,yTrafRot1); mi_setnodeprop('None',groupRot);
mi_selectnode(xArcRot1,yArcRot1); mi_setnodeprop('None',groupRot);
mi_clearselected;

mi_drawline(xArcRot2,yArcRot2,xTrafRot2,yTrafRot2);
mi_selectsegment(mean([xArcRot2 xTrafRot2]),mean([yArcRot2 yTrafRot2]));
mi_setsegmentprop('APg1', res_traf, 0, 0, groupRot);
mi_selectnode(xTrafRot2,yTrafRot2); mi_setnodeprop('None',groupRot);
mi_selectnode(xArcRot2,yArcRot2); mi_setnodeprop('None',groupRot);
mi_clearselected;

[xAirTrafRt,yAirTrafRt] = rot_point(r+1/6*g,0,0.5*ps*180/p*pi/180);
mi_addblocklabel(xAirTrafRt,yAirTrafRt);
mi_selectlabel(xAirTrafRt,yAirTrafRt);
mi_setblockprop('Air', 0, res_traf, 'None', 0, groupRot, 1);
mi_clearselected;

end