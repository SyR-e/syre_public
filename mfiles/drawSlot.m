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

function [geo,temp] = drawSlot(geo)

acs   = geo.acs;                        % slot opening [pu]
avv   = geo.win.avv;                    % widing matrix
lt    = geo.lt;                         % tooth length [mm]
wt    = geo.wt;                         % tooth width [mm]
p     = geo.p;                          % pole pairs number
q     = geo.q;                          % slot per pole per phase
R     = geo.R;                          % stator outer radius [mm]
r     = geo.r;                          % rotor outer radius [mm]
g     = geo.g;                          % airgap length [mm]
ttd   = geo.ttd;                        % tooth shoe thickness [mm]
tta   = geo.tta;                        % tooth shoe angle [°]
SFR   = geo.SFR;                        % fillet radius at the slot bottom [mm]
Qs    = geo.Qs;                         % number of simulated stator slots
n3ph  = geo.win.n3phase;                % number of 3phase sets
dis   = geo.win.liner;                  % liner thickness [mm]
pont0 = geo.pont0;                      % minimum mechanical tolerance [mm]
flagDivision = geo.statorYokeDivision;  % division line for end-winding extrusion (GalFer Challenge 2024)


slot_layer_pos = geo.win.slot_layer_pos;    % side-by-side flag winding

alpha_slot = 2*pi/(6*p*q*n3ph);  % angolo di passo cava
RSI        = r+g;           % r traferro statore
RSE        = R;             % r esterno statore
RS5        = r+g+lt;        % r esterno cave

slot_open_ang = acs*2*pi/(6*q*p*n3ph)/2;


% r eq for middle of slot computation
mr = tan(alpha_slot/2);
% x0 = RSI*cos(alpha_slot/2);
% y0 = RSI*sin(alpha_slot/2);

% if geo.parallel_slot
%  x2  = RSI+ttd;
%  y2  = x2*tan(slot_open_ang);
%  m23 = tan(pi/2-tta*pi/180);
%  q23 = y2-m23*x2;
%  x3 = (st/2-q23)/m23;
%  y3 = st/2;
%  %  [a23,b23,c23] = retta_mq2abc(m23,q23,nan);
% %  [x3,y3]= intersezione_tra_rette(0,0,y2,a23,b23,c23);
%  wt = abs(y3-mr*x3)/(1+mr^2)*2;
% end

% tooth side computation
% design like a line parallel to r, case of trapezoidal slot r2: y=m2x+q2
% explicit form
m2 = mr;
q2 = -wt/2*sqrt(1+mr^2);


[x1t,y1t]=intersezione_retta_circonferenza(0,0,RSI,m2,q2);


[x1,y1]=intersezione_retta_circonferenza(0,0,RSI,tan(slot_open_ang),0);
if (tan(y1./x1)>tan(y1t./x1t))
    x1=x1t;
    y1=y1t;
end

x2=x1+ttd;
y2=y1;

mtta=tan(pi/2-tta*pi/180);
qtta=y2-mtta*x2;
% ytta=mtta*xx+qtta;
[x3,y3]=intersezione_tra_rette(mtta,-1,qtta,m2,-1,q2);

st = y3*2;

% end of the slot
x6 = RSI+lt;
y6 = 0;
% LT2 position at the tooth
% q2 = y3;
% m2 = 0;
if geo.parallel_slot
    q2 = y3;
    m2 = 0;
end

[xLT2,yLT2]=intersezione_retta_circonferenza(0,0,(RSI+lt),m2,q2);
% bottom slot radius
[xRacSlot,yRacSlot,x5,y5,x4,y4]=cir_tg_2rette(x6,y6,xLT2,yLT2,xLT2,yLT2,x3,y3,SFR);
x4=x4(2);
y4=y4(2);
x5=x5(2);
y5=y5(2);
xRacSlot=xRacSlot(2);
yRacSlot=yRacSlot(2);

%% slot angle for SPM subdomain model
% geo.SlotAngle = atan(y4/x4)*2;                                              % bsa
% geo.SlotOpenAngle = slot_open_ang*2;                                        % boa

mm1 = (yLT2-y6)./(xLT2-x6);       % slope of line slot bottom
mm2 = (y3-yLT2)./(x3-xLT2);       % slope of line slot side
qq2 = yLT2-mm2*xLT2;
angle1 = atan(abs((mm2-mm1)./(1+mm1.*mm2)));        % angle between two lines
%% Chao limit fillet radius
if x5 > x6    % not beyond central of slot at slot bottom
    x5 = x6;
    y5 = y6;
    SFR = tan(angle1/2)*sqrt((yLT2-y5)^2+(xLT2-x5)^2);
    [xRacSlot,yRacSlot,x5,y5,x4,y4]=cir_tg_2rette(x6,y6,xLT2,yLT2,xLT2,yLT2,x3,y3,SFR);
    x4=x4(2);
    y4=y4(2);
    x5=x5(2);
    y5=y5(2);
    xRacSlot=xRacSlot(2);
    yRacSlot=yRacSlot(2);
end
if x4 < x2+(x6-x2)/2    % not beyond half of slot at slot side
    x4 = x2+(x6-x2)/2;
    y4 = m2*x4+q2;
    SFR = tan(angle1/2)*sqrt((yLT2-y4)^2+(xLT2-x4)^2);
    [xRacSlot,yRacSlot,x5,y5,x4,y4]=cir_tg_2rette(x6,y6,xLT2,yLT2,xLT2,yLT2,x3,y3,SFR);
    x4=x4(2);
    y4=y4(2);
    x5=x5(2);
    y5=y5(2);
    xRacSlot=xRacSlot(2);
    yRacSlot=yRacSlot(2);
end
area_corner = SFR^2 * (1./tan(angle1/2)-(pi-angle1)/2);           % redundant area at bottom slot area (made by slot side, slot bottom and circle)
% slot air (near air-gap)
xA1 = RSI;
yA1 = 0;
% slot air near copper
xA2 = x2;
yA2 = 0;
% external point
xE1 = RSE;
yE1 = 0;

xE2 = RSE*cos(alpha_slot/2);
yE2 = RSE*sin(alpha_slot/2);

% liner insertion (GalFer Challenge 2024, Repetto and Solimene)
x5_c=x5-dis;
y5_c=y5;

x4_c=x4+dis*sin(0.5*alpha_slot);
y4_c=y4-dis*cos(0.5*alpha_slot);
x2_c=x2+dis;
y2_c=y2;

xA2_c=xA2+dis;
yA2_c=0;

y6_c=0;
x6_c=x6-dis;

qtta_c=y2_c-mtta*x2_c;
m34=tan(alpha_slot/2);
q4=y4_c-tan(alpha_slot/2)*x4_c;
[x3_c,y3_c]=intersezione_tra_rette(mtta,-1,qtta_c,m34,-1,q4);

if y3_c<y2_c
    q4=y4_c-tan(alpha_slot/2)*x4_c;
    y3_c=tan(alpha_slot/2)*xA2_c+q4;
    x3_c=xA2_c;
    x2_c=x3_c;
    y2_c=y3_c;
end

%% slot area evaluation
if dis==0
    xTmp1 = x4-xRacSlot;
    yTmp1 = y4-yRacSlot;
    xTmp2 = x5-xRacSlot;
    yTmp2 = y5-yRacSlot;
    rRacc = abs(xTmp1+j*yTmp1);
    ang1  = angle(xTmp1+j*yTmp1);
    ang2  = angle(xTmp2+j*yTmp2);
    xRacc = xRacSlot+rRacc*cos(linspace(ang1,ang2,101));
    yRacc = yRacSlot+rRacc*sin(linspace(ang1,ang2,101));
    xArea = [xA2,x2,x3,xRacc,x6,xA2];
    yArea = [yA2,y2,y3,yRacc,y6,yA2];
    
else
    xTmp1 = x4_c-xRacSlot;
    yTmp1 = y4_c-yRacSlot;
    xTmp2 = x5_c-xRacSlot;
    yTmp2 = y5_c-yRacSlot;
    rRacc = abs(xTmp1+j*yTmp1);
    ang1  = angle(xTmp1+j*yTmp1);
    ang2  = angle(xTmp2+j*yTmp2);
    xRacc = xRacSlot+rRacc*cos(linspace(ang1,ang2,101));
    yRacc = yRacSlot+rRacc*sin(linspace(ang1,ang2,101));
    xArea = [xA2_c,x2_c,x3_c,xRacc,x6_c];
    yArea = [yA2_c,y2_c,y3_c,yRacc,y6_c];
end

%Aslot = 2*(polyarea(xArea,yArea)-area_corner);  % NB: the slot area should be computed from the slot matrix and not here...

Aslot = 2*(polyarea(xArea,yArea));

% Matrix describing lines and arches of half slot
materialCodes;

% codMatAirSta = 2;
% codMatFeSta  = 4;
% codMatCu  = 3;

% bottom side of the copper (airgap side)
if (xA2<xA1)
    CopperMat = [
        x3       0        x3 y3 NaN NaN  0 codMatCu 1
        x3       y3       x4 y4 NaN NaN  0 codMatCu 1
        xRacSlot yRacSlot x4 y4 x5  y5  -1 codMatCu 1
        x5       y5       x6 0  NaN NaN  0 codMatCu 1
        ];
    
    AirMat = [
        x1 y1 x2 y2 NaN NaN 0 codMatAirSta Qs+1
        x2 y2 x3 y3 NaN NaN 0 codMatAirSta Qs+1
        ];
    
else
    CopperMat = [
        x2       0        x2 y2 NaN NaN  0 codMatCu 1
        x2       y2       x3 y3 NaN NaN  0 codMatCu 1
        x3       y3       x4 y4 NaN NaN  0 codMatCu 1
        xRacSlot yRacSlot x4 y4 x5  y5  -1 codMatCu 1
        x5       y5       x6 0  NaN NaN  0 codMatCu 1
        ];
    
    AirMat = [
        x1 y1 x2 y2 NaN NaN 0 codMatAirSta Qs+1
        ];
end

if dis~=0
    CopperMat = [CopperMat;
        x2_c     0        x2_c y2_c NaN   NaN    0 codMatCu 2*Qs+1
        x2_c     y2_c     x3_c y3_c NaN   NaN    0 codMatCu 2*Qs+1
        x3_c     y3_c     x4_c y4_c NaN   NaN    0 codMatCu 2*Qs+1
        xRacSlot yRacSlot x4_c y4_c x5_c  y5_c  -1 codMatCu 2*Qs+1
        x5_c     y5_c     x6_c 0    NaN   NaN    0 codMatCu 2*Qs+1
        ];
end

% MIRROR of the half slot (adding the lines in the correct order)
CopperMatNeg = CopperMat;                           % create the second half slot
CopperMatNeg(:,[2,4,6]) = -CopperMat(:,[2,4,6]);    % mirror the points
CopperMatNeg = flipud(CopperMatNeg);                % invert the lines order (to close the object boundary)
indexArc = CopperMatNeg(:,7)~=0;                    % find the arc elements
CopperMatNeg(indexArc==0,1:4) = [CopperMatNeg(indexArc==0,3:4) CopperMatNeg(indexArc==0,1:2)];  % change the initial and final point of the lines
CopperMatNeg(indexArc==1,3:6) = [CopperMatNeg(indexArc==1,5:6) CopperMatNeg(indexArc==1,3:4)];  % change the initial and final point of the arcs

CopperMat = [CopperMat;CopperMatNeg];


AirMatNeg = AirMat;
AirMatNeg(:,[2,4,6]) = -AirMat(:,[2,4,6]);
AirMatNeg = flipud(AirMatNeg);
AirMatNeg(:,1:4) = [AirMatNeg(:,3:4) AirMatNeg(:,1:2)];

AirMat = [AirMat;AirMatNeg];

% winding layers and label centers

xCu = zeros(1,size(avv,1));
yCu = zeros(1,size(avv,1));

if (strcmp(slot_layer_pos,'side_by_side'))
    % if dis==0
        xCu(1) = (x2+x6)/2;
        yCu(1) = y5/2;
        xCu(2) = (x2+x6)/2;
        yCu(2) = -y5/2;
    % else
    %     xCu(1) = (x2_c+x6_c)/2;
    %     yCu(1) = y5_c/2;
    %     xCu(2) = (x2_c+x6_c)/2;
    %     yCu(2) = -y5_c/2;
    % end
    
    CopperMat = [CopperMat
        CopperMat(1,1)+dis CopperMat(1,2) x6-dis y6 NaN NaN 0 codMatCu 0
        ];
else
    % lines dividing the slot copper in 2 or more layers (bundles)
    Num_stat_strati=size(avv,1);
    if dis==0
        xini = CopperMat(1,1);
        xfin = xLT2;
        xstr = linspace(xini,xfin,Num_stat_strati+1);
        ystr = mm2*xstr+qq2;
    else
        xini = x2_c;
        xfin = x6_c;
        xstr = linspace(xini,xfin,Num_stat_strati+1);
        ystr = mm2*xstr+q4;
    end
    
    CopperMat = [CopperMat
        xstr(2:end-1)' ystr(2:end-1)' xstr(2:end-1)' -ystr(2:end-1)' NaN NaN 0 codMatCu 0
        ];
    
    xCu = xstr(1:end-1)+diff(xstr)/2;
    yCu = zeros(size(xCu));
end

if acs==0
    xAir = NaN;
    yAir = NaN;
    AirMat = [];
else
    xAir = (xA1+x2)/2;
    yAir = 0;
end

if dis~=0
    xAir = [xAir 0.5*(xA2+xA2_c)];
    yAir = [yAir 0];
end

xFe = [(RSE+RS5)/2 RS5+pont0/2];
yFe = [0 0];

if ~flagDivision
    xFe = xFe(1);
    yFe = yFe(1);
end

% xFe0 = RSE+pont0/2;
% yFe0 = 0;

slotMat = [CopperMat;AirMat];

% rotate in zero position
slotMat     = rotateMatrix(slotMat,alpha_slot/2);
[xCu,yCu]   = rot_point(xCu,yCu,alpha_slot/2);
[xAir,yAir] = rot_point(xAir,yAir,alpha_slot/2);
[xFe,yFe]   = rot_point(xFe,yFe,alpha_slot/2);
% [xFe0,yFe0] = rot_point(xFe0,yFe0,alpha_slot/2);

%% OUTPUT DATA
geo.Aslot = Aslot;
geo.SFR   = SFR;
% geo.wt = round(wt,2);
geo.st = round(st,2);

temp.xCu     = xCu;
temp.yCu     = yCu;
temp.xAir    = xAir;
temp.yAir    = yAir;
temp.xFe     = xFe;
temp.yFe     = yFe;
% temp.xFe0    = xFe0;
% temp.yFe0    = yFe0;
temp.CavaMat = slotMat;

% data for slot leakage computation
temp.d0 = geo.ttd;
temp.d1 = x3-x2;
temp.d2 = lt-temp.d0-temp.d1;
temp.c0 = 2*y2;
temp.c1 = 2*y3;
temp.c2 = 2*((xLT2-x6)^2+yLT2^2)^0.5;



