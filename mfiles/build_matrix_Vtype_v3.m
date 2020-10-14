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

function rotore = build_matrix_Vtype_v3(temp,geo)

nlay = geo.nlay;
PMdim = geo.PMdim(1,:);

% tangential ribs points
xpont      = temp.xpont;
ypont      = temp.ypont;
xxD1k      = temp.xxD1k;
yyD1k      = temp.yyD1k;
xxD2k      = temp.xxD2k;
yyD2k      = temp.yyD2k;
XcRibTraf1 = temp.XcRibTraf1;
YcRibTraf1 = temp.YcRibTraf1;
XcRibTraf2 = temp.XcRibTraf2;
YcRibTraf2 = temp.YcRibTraf2;
xTraf1     = temp.xTraf1;
yTraf1     = temp.yTraf1;
xTraf2     = temp.xTraf2;
yTraf2     = temp.yTraf2;

% radial ribs points
XpontRadSx    = temp.XpontRadSx;
YpontRadSx    = temp.YpontRadSx;
XpontRadDx    = temp.XpontRadDx;
YpontRadDx    = temp.YpontRadDx;
XpontRadBarDx = temp.XpontRadBarDx;
YpontRadBarDx = temp.YpontRadBarDx;
XpontRadBarSx = temp.XpontRadBarSx;
YpontRadBarSx = temp.YpontRadBarSx;


% magnet
xPMC1t = temp.xPMC1t;
yPMC1t = temp.yPMC1t;
xPMC2t = temp.xPMC2t;
yPMC2t = temp.yPMC2t;
xPMC1b = temp.xPMC1b;
yPMC1b = temp.yPMC1b;
xPMC2b = temp.xPMC2b;
yPMC2b = temp.yPMC2b;

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
    if (YpontRadSx(ii)~=0) % ponticello radiale
        rotore = [rotore
            XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadSx(ii)    YpontRadSx(ii)    NaN NaN 0 codMatAirRot indexEle
            XpontRadSx(ii)    YpontRadSx(ii)    XpontRadDx(ii)    YpontRadDx(ii)    NaN NaN 0 codMatAirRot indexEle
            XpontRadDx(ii)    YpontRadDx(ii)    XpontRadBarDx(ii) YpontRadBarDx(ii) NaN NaN 0 codMatAirRot indexEle
            ];
    else
        rotore = [rotore
            XpontRadBarSx(ii)   YpontRadBarSx(ii)   XpontRadBarDx(ii) YpontRadBarDx(ii) NaN               NaN               0 codMatAirRot indexEle
            ];
    end
    
    if isnan(xTraf1(ii)) % rounded end-barrier
        rotore = [rotore
            XpontRadBarDx(ii) YpontRadBarDx(ii) xxD2k(ii)         yyD2k(ii)         NaN       NaN       0 codMatAirRot indexEle
            XcRibTraf2(ii)    YcRibTraf2(ii)    xxD2k(ii)         yyD2k(ii)         xpont(ii) ypont(ii) 1 codMatAirRot indexEle
            XcRibTraf1(ii)    YcRibTraf1(ii)    xpont(ii)         ypont(ii)         xxD1k(ii) yyD1k(ii) 1 codMatAirRot indexEle
            xxD1k(ii)         yyD1k(ii)         XpontRadBarSx(ii) YpontRadBarSx(ii) NaN       NaN       0 codMatAirRot indexEle
            ];
    else
        rotore = [rotore
            XpontRadBarDx(ii) YpontRadBarDx(ii) xxD2k(ii)         yyD2k(ii)         NaN        NaN        0 codMatAirRot indexEle
            xxD2k(ii)         yyD2k(ii)         xTraf2(ii)        yTraf2(ii)        NaN        NaN        0 codMatAirRot indexEle
            0                 0                 xTraf2(ii)        yTraf2(ii)        xTraf1(ii) yTraf1(ii) 1 codMatAirRot indexEle
            xTraf1(ii)        yTraf1(ii)        xxD1k(ii)         yyD1k(ii)         NaN        NaN        0 codMatAirRot indexEle
            xxD1k(ii)         yyD1k(ii)         XpontRadBarSx(ii) YpontRadBarSx(ii) NaN        NaN        0 codMatAirRot indexEle
            ];
    end
    indexEle = indexEle+1;
    
    % Magnet
    if PMdim(ii)>0
        Mag = [Mag
            xPMC1b(ii) yPMC1b(ii) xPMC2b(ii) yPMC2b(ii) NaN NaN 0 codMatBar indexEle
            xPMC2b(ii) yPMC2b(ii) xPMC2t(ii) yPMC2t(ii) NaN NaN 0 codMatBar indexEle
            xPMC2t(ii) yPMC2t(ii) xPMC1t(ii) yPMC1t(ii) NaN NaN 0 codMatBar indexEle
            xPMC1t(ii) yPMC1t(ii) xPMC1b(ii) yPMC1b(ii) NaN NaN 0 codMatBar indexEle
            ];
        indexEle = indexEle+1;
    end
end

rotore=[rotore;Mag];



