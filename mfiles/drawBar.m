% Copyright 2022
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

function [geo,temp,mat] = drawBar(geo,mat,fem)

% This function draw a single rotor bar in zero position (as drawn in SyR-e
% for SyR machines), with labels matrix (BLKLABELSrot)

r         = geo.r;
wt        = geo.IM.wt;
lt        = geo.IM.lt;
acr       = geo.IM.acr;
ttd       = geo.IM.ttd;
Nbars     = geo.IM.Nbars;
filletTop = geo.IM.filletTop;
filletBot = geo.IM.filletBot;

alpha_slot = 2*pi/Nbars/2;
m = tan(alpha_slot);
q = 0;
% calcolo mezza cava
[xOut,yOut] = intersezione_retta_circonferenza(0,0,r,m,q);
[aOut,bOut,cOut] = retta_per_2pti(0,0,xOut,yOut);
[a,b,c] = calc_retta_offset(aOut,bOut,cOut,wt/2);

m1 = tan(alpha_slot*acr);
q1 = 0;
% lati apertura cava
[x1,y1] = intersezione_retta_circonferenza(0,0,r,m1,q1);
[x2,y2] = intersezione_retta_circonferenza(0,0,r-ttd,m1,q1);


% calcolo raccordo superiore
[ap,bp,cp] = calc_retta_offset(a,b,c,filletTop);
[xTmp,yTmp] = calc_int_retta_circ_gen(0,0,r-ttd,ap,bp,cp);
x31 = xTmp(1);
y31 = yTmp(1);
[xTmp,yTmp] = calc_int_retta_circ_gen(0,0,r-ttd-filletTop,a,b,c);
x32 = xTmp(1);
y32 = yTmp(1);
[xTmp,yTmp] = calc_int_retta_circ_gen(0,0,r-ttd-filletTop,ap,bp,cp);
x30 = xTmp(1);
y30 = yTmp(1);
% calcolo fondo cava
[ap,bp,cp] = calc_retta_offset(a,b,c,filletBot);
aB = -1;
bB = 0;
cB = r-lt;
aBp = aB;
bBp = bB;
cBp = r-lt+filletBot;
[x40,y40] = intersezione_tra_rette(ap,bp,cp,aBp,bBp,cBp);
[x41,y41] = intersezione_tra_rette(a,b,c,aBp,bBp,cBp);
[x42,y42] = intersezione_tra_rette(ap,bp,cp,aB,bB,cB);


% creazione della matrice rotore
materialCodes;
rotor = [
        x1   y1  x2   y2  NaN  NaN 0 codMatCuRot 1
        0    0   x2   y2  x31  y31 1 codMatCuRot 1
        x30  y30 x31  y31 x32  y32 1 codMatCuRot 1
        x32  y32 x41  y41 NaN  NaN 0 codMatCuRot 1
        x40  y40 x41  y41 x42  y42 1 codMatCuRot 1
        x42  y42 x42  0   NaN  NaN 0 codMatCuRot 1
        x42  0   x42 -y42 NaN  NaN 0 codMatCuRot 1
        x40 -y40 x42 -y42 x41 -y41 1 codMatCuRot 1
        x41 -y41 x32 -y32 NaN  NaN 0 codMatCuRot 1
        x30 -y30 x32 -y32 x31 -y31 1 codMatCuRot 1
        0    0   x31 -y31 x2  -y2  1 codMatCuRot 1
        x2  -y2  x1  -y1  NaN  NaN 0 codMatCuRot 1
    ];

% definizione dei labels
xy = [r-ttd-lt/2 0 codMatCuRot fem.res 1 0 0 0];

% rotazione in posizione 0
rotor = rotateMatrix(rotor,pi/Nbars);
rotor = checkPlotMatrix(rotor,1e-9);

[xtemp,ytemp] = rot_point(xy(:,1),xy(:,2),pi/Nbars);
xy = [xtemp ytemp xy(:,3:end)];

% output
temp.x1  = x1;
temp.y1  = y1;
temp.x2  = x2;
temp.y2  = y2;
temp.x30 = x30;
temp.y30 = y30;
temp.x31 = x31;
temp.y31 = y31;
temp.x32 = x32;
temp.y32 = y32;
temp.x40 = x40;
temp.y40 = y40;
temp.x41 = x41;
temp.y41 = y41;
temp.x42 = x42;
temp.y42 = y42;

geo.rotor = rotor;
geo.BLKLABELS.rotore.xy = xy;




















