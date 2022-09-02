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

function rotore = build_matrix_Circ_dx(temp,geo)

x0    = geo.x0;
PMdim = geo.PMdim(1,:);
nlay  = geo.nlay;

RotorFilletTan1 = geo.RotorFilletTan1;
RotorFilletTan2 = geo.RotorFilletTan2;

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

XpontRadSx    = temp.XpontRadSx;
YpontRadSx    = temp.YpontRadSx;
XpontRadDx    = temp.XpontRadDx;
YpontRadDx    = temp.YpontRadDx;
XpontRadBarDx = temp.XpontRadBarDx;
YpontRadBarDx = temp.YpontRadBarDx;
XpontRadBarSx = temp.XpontRadBarSx;
YpontRadBarSx = temp.YpontRadBarSx;

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

if isfield(temp,'xcRac1')
    xcRac1 = temp.xcRac1;
    ycRac1 = temp.ycRac1;
    xxE1k =temp.xxE1k;
    yyE1k =temp.yyE1k;
    xcRac2 = temp.xcRac2;
    ycRac2 = temp.ycRac2;
    xxE2k =temp.xxE2k;
    yyE2k =temp.yyE2k;
end

xS01k = temp.xS01k;
yS01k = temp.yS01k;
xS02k = temp.xS02k;
yS02k = temp.yS02k;


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
    if (YpontRadSx(ii)~=0)
%         rotore = [rotore
%             XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadSx(ii)    YpontRadSx(ii)    NaN NaN 0 codMatAirRot indexEle
%             XpontRadSx(ii)    YpontRadSx(ii)    XpontRadDx(ii)    YpontRadDx(ii)    NaN NaN 0 codMatAirRot indexEle
%             XpontRadDx(ii)    YpontRadDx(ii)    XpontRadBarDx(ii) YpontRadBarDx(ii) NaN NaN 0 codMatAirRot indexEle
%             ];
        rotore = [rotore
            xS01k(ii)         yS01k(ii)      XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadSx(ii)    YpontRadSx(ii)    1 codMatAirRot indexEle
            XpontRadSx(ii)    YpontRadSx(ii) XpontRadDx(ii)    YpontRadDx(ii)    NaN               NaN               0 codMatAirRot indexEle
            xS02k(ii)         yS02k(ii)      XpontRadDx(ii)    YpontRadDx(ii)    XpontRadBarDx(ii) YpontRadBarDx(ii) 1 codMatAirRot indexEle
            ];
    else
        rotore = [rotore
            XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadBarDx(ii) YpontRadBarDx(ii) NaN NaN 0 codMatAirRot indexEle
            ];
    end
    
    if ~(isfinite(RotorFilletTan1(ii)) && isfinite(RotorFilletTan2(ii)))
        rotore = [rotore
            x0             0              XpontRadBarDx(ii) YpontRadBarDx(ii) xxD2k(ii)         yyD2k(ii)         -1 codMatAirRot indexEle
            XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii)         yyD2k(ii)         xpont(ii)         ypont(ii)         +1 codMatAirRot indexEle
            XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii)         ypont(ii)         xxD1k(ii)         yyD1k(ii)         +1 codMatAirRot indexEle
            x0             0              xxD1k(ii)         yyD1k(ii)         XpontRadBarSx(ii) YpontRadBarSx(ii) +1 codMatAirRot indexEle
            ];

        indexEle = indexEle+1;

    else

        % Rotor fillet
        rotore = [rotore
            x0             0              XpontRadBarDx(ii) YpontRadBarDx(ii) xxD2k(ii)         yyD2k(ii)         -1 codMatAirRot indexEle
            xcRac2(ii)     ycRac2(ii)     xxD2k(ii)         yyD2k(ii)         xxE2k(ii)         yyE2k(ii)         +1 codMatAirRot indexEle
            0              0              xxE2k(ii)         yyE2k(ii)         xxE1k(ii)         yyE1k(ii)         +1 codMatAirRot indexEle
            xcRac1(ii)     ycRac1(ii)     xxE1k(ii)         yyE1k(ii)         xxD1k(ii)         yyD1k(ii)         +1 codMatAirRot indexEle
            x0             0              xxD1k(ii)         yyD1k(ii)         XpontRadBarSx(ii) YpontRadBarSx(ii) +1 codMatAirRot indexEle

            ];

        indexEle = indexEle+1;
    end

    % Magnets
    if PMdim(ii)>0
        Mag = [Mag
            xPMC1b(ii) yPMC1b(ii) xPMC2b(ii) yPMC2b(ii) NaN        NaN         0 codMatBar indexEle
            x0         0          xPMC2b(ii) yPMC2b(ii) xPMC2t(ii) yPMC2t(ii) -1 codMatBar indexEle
            xPMC2t(ii) yPMC2t(ii) xPMC1t(ii) yPMC1t(ii) NaN        NaN         0 codMatBar indexEle
            x0         0          xPMC1t(ii) yPMC1t(ii) xPMC1b(ii) yPMC1b(ii) +1 codMatBar indexEle
            ];
        
        indexEle = indexEle+1;
        
        Mag = [Mag
            xPME1b(ii) yPME1b(ii) xPME2b(ii) yPME2b(ii) NaN        NaN         0 codMatBar indexEle
            x0         0          xPME2b(ii) yPME2b(ii) xPME2t(ii) yPME2t(ii) -1 codMatBar indexEle
            xPME2t(ii) yPME2t(ii) xPME1t(ii) yPME1t(ii) NaN        NaN         0 codMatBar indexEle
            x0         0          xPME1t(ii) yPME1t(ii) xPME1b(ii) yPME1b(ii) +1 codMatBar indexEle
            ];
        
        indexEle = indexEle+1;
    end
end

rotore = [rotore;Mag];
