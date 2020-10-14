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

function rotore = build_matrix_ISeg(temp,geo)

PMdimC   = geo.PMdim(1,:);
PMdimE   = geo.PMdim(2,:);
nlay     = geo.nlay;
pontR    = geo.pontR;
splitRib = geo.radial_ribs_split;

B1k   = temp.B1k;
B2k   = temp.B2k;
xpont = temp.xpont;
ypont = temp.ypont;

XcRibTraf1 = temp.XcRibTraf1;
YcRibTraf1 = temp.YcRibTraf1;
XcRibTraf2 = temp.XcRibTraf2;
YcRibTraf2 = temp.YcRibTraf2;

XpBar1 = temp.XpBar1;
YpBar1 = temp.YpBar1;
XpBar2 = temp.XpBar2;
YpBar2 = temp.YpBar2;

xxD1k = temp.xxD1k;
yyD1k = temp.yyD1k;
xxD2k = temp.xxD2k;
yyD2k = temp.yyD2k;

XpontRadSx    = temp.XpontRadSx;
YpontRadSx    = temp.YpontRadSx;
XpontRadDx    = temp.XpontRadDx;
YpontRadDx    = temp.YpontRadDx;
XpontRadBarDx = temp.XpontRadBarDx;
YpontRadBarDx = temp.YpontRadBarDx;
XpontRadBarSx = temp.XpontRadBarSx;
YpontRadBarSx = temp.YpontRadBarSx;

XpontSplitBarSx = temp.XpontSplitBarSx;
YpontSplitBarSx = temp.YpontSplitBarSx;
XpontSplitBarDx = temp.XpontSplitBarDx;
YpontSplitBarDx = temp.YpontSplitBarDx;
XpontSplitDx    = temp.XpontSplitDx;
YpontSplitDx    = temp.YpontSplitDx;
XpontSplitSx    = temp.XpontSplitSx;
YpontSplitSx    = temp.YpontSplitSx;

xPMC1b = temp.xPMC1b;
yPMC1b = temp.yPMC1b;
xPMC1t = temp.xPMC1t;
yPMC1t = temp.yPMC1t;
xPMC2b = temp.xPMC2b;
yPMC2b = temp.yPMC2b;
xPMC2t = temp.xPMC2t;
yPMC2t = temp.yPMC2t;
xPME1b = temp.xPME1b;
yPME1b = temp.yPME1b;
xPME1t = temp.xPME1t;
yPME1t = temp.yPME1t;
xPME2b = temp.xPME2b;
yPME2b = temp.yPME2b;
xPME2t = temp.xPME2t;
yPME2t = temp.yPME2t;


% This function build the rotor matrix (defining the geometry). Each half
% barrier and PM is composed clockwise, from the bottom-left corner.
% Each line of the matrix represent a line or an arch, according to the
% following standard:
% - line:
%        [x1 y1 x2 y2 NaN NaN 0 codMat numEle]
% - arc:
%        [x0 y0 x1 y1  x2  y2 direction codMat numEle]
% where:
% - (x0,y0)   --> center of the arc
% - (x1,y1)   --> initial point of the line/arc
% - (x2,y2)   --> final point of the line/arc
% - direction -->  0 --> line
%                 +1 --> arc, from 1 to 2, clockwise
%                 -1 --> arc, from 2 to 1, counter-clockwise
% codMat      --> 1 --> air
%                 6 --> PM
% numEle      --> element number of holes in the rotor lamination or PM
%                 (two different lists)

materialCodes;
% codMatAirRot = 1;
% codMatBar  = 6;

rotore  = [];
Mag     = [];

indexEle = 1;

for ii=1:nlay
    % rotor matrix
    if pontR(ii)~=0
        if splitRib
            rotore = [rotore
                XpontRadBarSx(ii)     YpontRadBarSx(ii)     XpontRadBarDx(ii)     YpontRadBarDx(ii)     NaN NaN 0 codMatAirRot indexEle
                XpontRadBarDx(ii)     YpontRadBarDx(ii)     XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) NaN NaN 0 codMatAirRot indexEle
                XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) XpontSplitDx(2,ii)    YpontSplitDx(2,ii)    NaN NaN 0 codMatAirRot indexEle
                XpontSplitDx(2,ii)    YpontSplitDx(2,ii)    XpontSplitSx(2,ii)    YpontSplitSx(2,ii)    NaN NaN 0 codMatAirRot indexEle
                XpontSplitSx(2,ii)    YpontSplitSx(2,ii)    XpontSplitBarSx(2,ii) YpontSplitBarSx(2,ii) NaN NaN 0 codMatAirRot indexEle
                XpontSplitBarSx(2,ii) YpontSplitBarSx(2,ii) XpontRadBarSx(ii)     YpontRadBarSx(ii)     NaN NaN 0 codMatAirRot indexEle
                ];
            indexEle = indexEle+1;
            
            rotore = [rotore
                XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    NaN       NaN       0 codMatAirRot indexEle
                XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    NaN       NaN       0 codMatAirRot indexEle
                XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) NaN       NaN       0 codMatAirRot indexEle
                XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) xxD2k(ii)             yyD2k(ii)             NaN       NaN       0 codMatAirRot indexEle
                XcRibTraf2(ii)        YcRibTraf2(ii)        xxD2k(ii)             yyD2k(ii)             xpont(ii) ypont(ii) 1 codMatAirRot indexEle
                XcRibTraf1(ii)        YcRibTraf1(ii)        xpont(ii)             ypont(ii)             xxD1k(ii) yyD1k(ii) 1 codMatAirRot indexEle
                xxD1k(ii)             yyD1k(ii)             XpBar1(ii)            YpBar1(ii)            NaN       NaN       0 codMatAirRot indexEle
                XpBar1(ii)        YpBar1(ii)                XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) NaN       NaN       0 codMatAirRot indexEle
                ];
            indexEle = indexEle+1;
        else
            rotore = [rotore
                XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadSx(ii)    YpontRadSx(ii)    NaN       NaN       0 codMatAirRot indexEle
                XpontRadSx(ii)    YpontRadSx(ii)    XpontRadDx(ii)    YpontRadDx(ii)    NaN       NaN       0 codMatAirRot indexEle
                XpontRadDx(ii)    YpontRadDx(ii)    XpontRadBarDx(ii) YpontRadBarDx(ii) NaN       NaN       0 codMatAirRot indexEle
                XpontRadBarDx(ii) YpontRadBarDx(ii) XpBar2(ii)        YpBar2(ii)        NaN       NaN       0 codMatAirRot indexEle
                XpBar2(ii)        YpBar2(ii)        xxD2k(ii)         yyD2k(ii)         NaN       NaN       0 codMatAirRot indexEle
                XcRibTraf2(ii)    YcRibTraf2(ii)    xxD2k(ii)         yyD2k(ii)         xpont(ii) ypont(ii) 1 codMatAirRot indexEle
                XcRibTraf1(ii)    YcRibTraf1(ii)    xpont(ii)         ypont(ii)         xxD1k(ii) yyD1k(ii) 1 codMatAirRot indexEle
                xxD1k(ii)         yyD1k(ii)         XpBar1(ii)        YpBar1(ii)        NaN       NaN       0 codMatAirRot indexEle
                XpBar1(ii)        YpBar1(ii)        XpontRadBarSx(ii) YpontRadBarSx(ii) NaN       NaN       0 codMatAirRot indexEle
                ];
            indexEle = indexEle+1;
        end
    else
        rotore = [rotore
            XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadBarDx(ii) YpontRadBarDx(ii) NaN       NaN       0 codMatAirRot indexEle
            XpontRadBarDx(ii) YpontRadBarDx(ii) XpBar2(ii)        YpBar2(ii)        NaN       NaN       0 codMatAirRot indexEle
            XpBar2(ii)        YpBar2(ii)        xxD2k(ii)         yyD2k(ii)         NaN       NaN       0 codMatAirRot indexEle
            XcRibTraf2(ii)    YcRibTraf2(ii)    xxD2k(ii)         yyD2k(ii)         xpont(ii) ypont(ii) 1 codMatAirRot indexEle
            XcRibTraf1(ii)    YcRibTraf1(ii)    xpont(ii)         ypont(ii)         xxD1k(ii) yyD1k(ii) 1 codMatAirRot indexEle
            xxD1k(ii)         yyD1k(ii)         XpBar1(ii)        YpBar1(ii)        NaN       NaN       0 codMatAirRot indexEle
            XpBar1(ii)        YpBar1(ii)        XpontRadBarSx(ii) YpontRadBarSx(ii) NaN       NaN       0 codMatAirRot indexEle
            ];
        
        indexEle = indexEle+1;
    end
    
    % PM
    if PMdimC(ii)>0
        Mag = [Mag
            xPMC1b(ii) yPMC1b(ii) xPMC2b(ii) yPMC2b(ii) NaN NaN 0 codMatBar indexEle
            xPMC2b(ii) yPMC2b(ii) xPMC2t(ii) yPMC2t(ii) NaN NaN 0 codMatBar indexEle
            xPMC2t(ii) yPMC2t(ii) xPMC1t(ii) yPMC1t(ii) NaN NaN 0 codMatBar indexEle
            xPMC1t(ii) yPMC1t(ii) xPMC1b(ii) yPMC1b(ii) NaN NaN 0 codMatBar indexEle
            ];
        indexEle = indexEle+1;
    end
    
    if PMdimE(ii)>0
        Mag = [Mag
            xPME1b(ii) yPME1b(ii) xPME2b(ii) yPME2b(ii) NaN NaN 0 codMatBar indexEle
            xPME2b(ii) yPME2b(ii) xPME2t(ii) yPME2t(ii) NaN NaN 0 codMatBar indexEle
            xPME2t(ii) yPME2t(ii) xPME1t(ii) yPME1t(ii) NaN NaN 0 codMatBar indexEle
            xPME1t(ii) yPME1t(ii) xPME1b(ii) yPME1b(ii) NaN NaN 0 codMatBar indexEle
            ];
        indexEle = indexEle+1;
    end
    
end

% Inserimento matrice Mag
rotore=[rotore;Mag];



%% OLD FUNCTION
Bx0=temp.Bx0;
B1k=temp.B1k;
B2k=temp.B2k;
xpont=temp.xpont;
ypont=temp.ypont;
xxD1k=temp.xxD1k;
yyD1k=temp.yyD1k;
xxD2k=temp.xxD2k;
yyD2k=temp.yyD2k;
XpBar1=temp.XpBar1;
YpBar1=temp.YpBar1;
XpBar2=temp.XpBar2;
YpBar2=temp.YpBar2;

xTraf1=temp.xTraf1;
xTraf2=temp.xTraf2;
yTraf1=temp.yTraf1;
yTraf2=temp.yTraf2;

arcLayTraf1=temp.arcLayTraf1;
arcLayTraf2=temp.arcLayTraf2;

XcRibTraf1=temp.XcRibTraf1;
YcRibTraf1=temp.YcRibTraf1;
XcRibTraf2=temp.XcRibTraf2;
YcRibTraf2=temp.YcRibTraf2;
% radial ribs coordinate
XpontRadSx=temp.XpontRadSx;
YpontRadSx=temp.YpontRadSx;
XpontRadDx=temp.XpontRadDx;
YpontRadDx=temp.YpontRadDx;
XpontRadBarDx=temp.XpontRadBarDx;
YpontRadBarDx=temp.YpontRadBarDx;
XpontRadBarSx=temp.XpontRadBarSx;
YpontRadBarSx=temp.YpontRadBarSx;
% Additional point for magnet
XpMag1B1=temp.XpMag1B1;
YpMag1B1=temp.YpMag1B1;

% Error mex no linear barrier boundary
error_mex=temp.error_mex;
% elementi per completamento punti
Mag=temp.Mag;

xob1pt1=temp.xob1pt1;
yob1pt1=temp.yob1pt1;
xob1pt2=temp.xob1pt2;
yob1pt2=temp.yob1pt2;
xob2pt1=temp.xob2pt1;
yob2pt1=temp.yob2pt1;
xob2pt2=temp.xob2pt2;
yob2pt2=temp.yob2pt2;


rotore=[];

for ii=1:geo.nlay
    
    if (YpontRadSx(ii)~=0)
        rotore=[rotore;XpontRadSx(ii),YpontRadSx(ii),XpontRadDx(ii),YpontRadDx(ii),NaN,NaN,0;...
            XpontRadSx(ii),YpontRadSx(ii),XpontRadBarSx(ii) YpontRadBarSx(ii),NaN,NaN,0;...
            XpontRadDx(ii),YpontRadDx(ii),XpontRadBarDx(ii) YpontRadBarDx(ii),NaN,NaN,0];
    else
        rotore=[rotore;
            B1k(ii),0,B2k(ii),0,NaN,NaN,0];
    end
    if (error_mex(ii)==0)
        if ii==1
            rotore=[rotore;
                XpontRadBarSx(ii) YpontRadBarSx(ii) XpBar1(ii) YpBar1(ii) NaN NaN 0;
                XpontRadBarDx(ii) YpontRadBarDx(ii) XpBar2(ii) YpBar2(ii) NaN NaN 0;
                XpBar1(ii) YpBar1(ii) xxD1k(ii) yyD1k(ii) NaN NaN 0;
                XpBar2(ii) YpBar2(ii) xxD2k(ii) yyD2k(ii) NaN NaN 0];
            rotore=[rotore;
                XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii) ypont(ii) xxD1k(ii) yyD1k(ii) 1];
            rotore=[rotore;
                XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii) yyD2k(ii) xpont(ii) ypont(ii) 1];
        else
            rotore=[rotore;
                XpontRadBarSx(ii) YpontRadBarSx(ii) XpBar1(ii) YpBar1(ii) NaN NaN 0;
                XpontRadBarDx(ii) YpontRadBarDx(ii) XpBar2(ii) YpBar2(ii) NaN NaN 0;
                XpBar1(ii) YpBar1(ii) xob1pt1(ii) yob1pt1(ii) NaN NaN 0;
                xob1pt1(ii) yob1pt1(ii) xob1pt2(ii) yob1pt2(ii) NaN NaN 0;
                xob1pt2(ii) yob1pt2(ii) xxD1k(ii) yyD1k(ii) NaN NaN 0;
                XpBar2(ii) YpBar2(ii) xob2pt2(ii) yob2pt2(ii) NaN NaN 0;
                xob2pt2(ii) yob2pt2(ii) xxD2k(ii) yyD2k(ii) NaN NaN 0];
            rotore=[rotore;
                XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii) ypont(ii) xxD1k(ii) yyD1k(ii) 1];
            rotore=[rotore;
                XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii) yyD2k(ii) xpont(ii) ypont(ii) 1];
            
        end
        
    else
        rotore=[rotore;
            XcRibTraf1(ii) YcRibTraf1(ii) xxD2k(ii) yyD2k(ii) xxD1k(ii) yyD1k(ii) 1];
        rotore=[rotore;
            XcRibTraf2(ii) YcRibTraf2(ii) xxD1k(ii) yyD1k(ii) xxD2k(ii) yyD2k(ii) 1];
        
    end
    
    
    
end

%% Inserimento matrice Mag

rotore=[rotore;Mag];

end

