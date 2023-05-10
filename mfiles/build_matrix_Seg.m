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

function rotore = build_matrix_Seg(temp,geo)

PMdimC   = geo.PMdim(1,:);
PMdimE   = geo.PMdim(2,:);
PMNcE    = geo.PMNc(2,:);
PMNcC    = geo.PMNc(1,:);
nlay     = geo.nlay;
pontR    = geo.pontR;

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

xS01k = temp.xS01k;
yS01k = temp.yS01k;
xS02k = temp.xS02k;
yS02k = temp.yS02k;
xI01k = temp.xI01k;
yI01k = temp.yI01k;
xI02k = temp.xI02k;
yI02k = temp.yI02k;

flag_segV      = temp.flag_segV;
flagPosSplitDx = temp.flagPosSplitDx;
flagPosSplitSx = temp.flagPosSplitSx;

radial_ribs_split = geo.radial_ribs_split;

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
yPME2t = temp.yPME2t;
xPME2s = temp.xPME2s;
yPME2s = temp.yPME2s;
xPME1s = temp.xPME1s;
yPME1s = temp.yPME1s;
xPMC2s = temp.xPMC2s;
yPMC2s = temp.yPMC2s;
xPMC1s = temp.xPMC1s;
yPMC1s = temp.yPMC1s;

xC1k  = temp.xC1k;
yC1k  = temp.yC1k;
xC2k  = temp.xC2k;
yC2k  = temp.yC2k;
xC3k  = temp.xC3k;
yC3k  = temp.yC3k;
xC4k  = temp.xC4k;
yC4k  = temp.yC4k;
xC01k = temp.xC01k;
yC01k = temp.yC01k;
xC02k = temp.xC02k;
yC02k = temp.yC02k;

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
    % first check: radial ribs?
    if pontR(ii)==0
        % tangential ribs shape
        if isnan(xC3k(ii)) % classical tangential rib
            rotore = [rotore
                XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadBarDx(ii) YpontRadBarDx(ii) NaN       NaN       0 codMatAirRot indexEle
                XpontRadBarDx(ii) YpontRadBarDx(ii) XpBar2(ii)        YpBar2(ii)        NaN       NaN       0 codMatAirRot indexEle
                XpBar2(ii)        YpBar2(ii)        xxD2k(ii)         yyD2k(ii)         NaN       NaN       0 codMatAirRot indexEle
                XcRibTraf2(ii)    YcRibTraf2(ii)    xxD2k(ii)         yyD2k(ii)         xpont(ii) ypont(ii) 1 codMatAirRot indexEle
                XcRibTraf1(ii)    YcRibTraf1(ii)    xpont(ii)         ypont(ii)         xxD1k(ii) yyD1k(ii) 1 codMatAirRot indexEle
                xxD1k(ii)         yyD1k(ii)         XpBar1(ii)        YpBar1(ii)        NaN       NaN       0 codMatAirRot indexEle
                XpBar1(ii)        YpBar1(ii)        XpontRadBarSx(ii) YpontRadBarSx(ii) NaN       NaN       0 codMatAirRot indexEle
                ];
        else % straight tangential rib
            rotore = [rotore
                XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadBarDx(ii) YpontRadBarDx(ii) NaN       NaN       0 codMatAirRot indexEle
                XpontRadBarDx(ii) YpontRadBarDx(ii) XpBar2(ii)        YpBar2(ii)        NaN       NaN       0 codMatAirRot indexEle
                XpBar2(ii)        YpBar2(ii)        xxD2k(ii)         yyD2k(ii)         NaN       NaN       0 codMatAirRot indexEle
                xC02k(ii)         yC02k(ii)         xxD2k(ii)         yyD2k(ii)         xC3k(ii)  yC3k(ii)  1 codMatAirRot indexEle
                0                 0                 xC3k(ii)          yC3k(ii)          xC2k(ii)  yC2k(ii)  1 codMatAirRot indexEle
                xC01k(ii)         yC01k(ii)         xC2k(ii)          yC2k(ii)          xxD1k(ii) yyD1k(ii) 1 codMatAirRot indexEle
                xxD1k(ii)         yyD1k(ii)         XpBar1(ii)        YpBar1(ii)        NaN       NaN       0 codMatAirRot indexEle
                XpBar1(ii)        YpBar1(ii)        XpontRadBarSx(ii) YpontRadBarSx(ii) NaN       NaN       0 codMatAirRot indexEle
                ];
        end
        indexEle = indexEle+1;
    else
        if radial_ribs_split(ii)==0 % radial ribs single
            if isnan(xC3k(ii)) % classical tangential rib
                rotore = [rotore
                    xS01k(ii)         yS01k(ii)         XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadSx(ii)    YpontRadSx(ii)    1 codMatAirRot indexEle
                    XpontRadSx(ii)    YpontRadSx(ii)    XpontRadDx(ii)    YpontRadDx(ii)    NaN               NaN               0 codMatAirRot indexEle
                    xS02k(ii)         yS02k(ii)         XpontRadDx(ii)    YpontRadDx(ii)    XpontRadBarDx(ii) YpontRadBarDx(ii) 1 codMatAirRot indexEle
                    XpontRadBarDx(ii) YpontRadBarDx(ii) XpBar2(ii)        YpBar2(ii)        NaN               NaN               0 codMatAirRot indexEle
                    XpBar2(ii)        YpBar2(ii)        xxD2k(ii)         yyD2k(ii)         NaN               NaN               0 codMatAirRot indexEle
                    XcRibTraf2(ii)    YcRibTraf2(ii)    xxD2k(ii)         yyD2k(ii)         xpont(ii)         ypont(ii)         1 codMatAirRot indexEle
                    XcRibTraf1(ii)    YcRibTraf1(ii)    xpont(ii)         ypont(ii)         xxD1k(ii)         yyD1k(ii)         1 codMatAirRot indexEle
                    xxD1k(ii)         yyD1k(ii)         XpBar1(ii)        YpBar1(ii)        NaN               NaN               0 codMatAirRot indexEle
                    XpBar1(ii)        YpBar1(ii)        XpontRadBarSx(ii) YpontRadBarSx(ii) NaN               NaN               0 codMatAirRot indexEle
                    ];
            else % straight tangential rib
                rotore = [rotore
                    xS01k(ii)         yS01k(ii)         XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadSx(ii)    YpontRadSx(ii)    1 codMatAirRot indexEle
                    XpontRadSx(ii)    YpontRadSx(ii)    XpontRadDx(ii)    YpontRadDx(ii)    NaN               NaN               0 codMatAirRot indexEle
                    xS02k(ii)         yS02k(ii)         XpontRadDx(ii)    YpontRadDx(ii)    XpontRadBarDx(ii) YpontRadBarDx(ii) 1 codMatAirRot indexEle
                    XpontRadBarDx(ii) YpontRadBarDx(ii) XpBar2(ii)        YpBar2(ii)        NaN               NaN               0 codMatAirRot indexEle
                    XpBar2(ii)        YpBar2(ii)        xxD2k(ii)         yyD2k(ii)         NaN               NaN               0 codMatAirRot indexEle
                    xC02k(ii)         yC02k(ii)         xxD2k(ii)         yyD2k(ii)         xC3k(ii)          yC3k(ii)          1 codMatAirRot indexEle
                    0                 0                 xC3k(ii)          yC3k(ii)          xC2k(ii)          yC2k(ii)          1 codMatAirRot indexEle
                    xC01k(ii)         yC01k(ii)         xC2k(ii)          yC2k(ii)          xxD1k(ii)         yyD1k(ii)         1 codMatAirRot indexEle
                    xxD1k(ii)         yyD1k(ii)         XpBar1(ii)        YpBar1(ii)        NaN               NaN               0 codMatAirRot indexEle
                    XpBar1(ii)        YpBar1(ii)        XpontRadBarSx(ii) YpontRadBarSx(ii) NaN               NaN               0 codMatAirRot indexEle
                    ];
            end
            indexEle = indexEle+1;
        else
            % central section
            if flagPosSplitSx(2,ii)==0 % bottom side of the ribs on the central section
                rotore = [rotore
                    XpontRadBarSx(ii)     YpontRadBarSx(ii)     XpontRadBarDx(ii)     YpontRadBarDx(ii)     NaN                   NaN                   0 codMatAirRot indexEle
                    XpontRadBarDx(ii)     YpontRadBarDx(ii)     XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) NaN                   NaN                   0 codMatAirRot indexEle
                    xI02k(ii)             yI02k(ii)             XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) XpontSplitDx(2,ii)    YpontSplitDx(2,ii)    1 codMatAirRot indexEle
                    XpontSplitDx(2,ii)    YpontSplitDx(2,ii)    XpontSplitSx(2,ii)    YpontSplitSx(2,ii)    NaN                   NaN                   0 codMatAirRot indexEle
                    xI01k(ii)             yI01k(ii)             XpontSplitSx(2,ii)    YpontSplitSx(2,ii)    XpontSplitBarSx(2,ii) YpontSplitBarSx(2,ii) 1 codMatAirRot indexEle
                    XpontSplitBarSx(2,ii) YpontSplitBarSx(2,ii) XpontRadBarSx(ii)     YpontRadBarSx(ii)     NaN                   NaN                   0 codMatAirRot indexEle
                    ];
                indexEle = indexEle+1;
            else
                rotore = [rotore
                    XpontRadBarSx(ii)     YpontRadBarSx(ii)     XpontRadBarDx(ii)     YpontRadBarDx(ii)     NaN                   NaN                   0 codMatAirRot indexEle
                    XpontRadBarDx(ii)     YpontRadBarDx(ii)     XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) NaN                   NaN                   0 codMatAirRot indexEle
                    xI02k(ii)             yI02k(ii)             XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) XpontSplitDx(2,ii)    YpontSplitDx(2,ii)    1 codMatAirRot indexEle
                    XpontSplitDx(2,ii)    YpontSplitDx(2,ii)    XpontSplitSx(2,ii)    YpontSplitSx(2,ii)    NaN                   NaN                   0 codMatAirRot indexEle
                    xI01k(ii)             yI01k(ii)             XpontSplitSx(2,ii)    YpontSplitSx(2,ii)    XpontSplitBarSx(2,ii) YpontSplitBarSx(2,ii) 1 codMatAirRot indexEle
                    XpontSplitBarSx(2,ii) YpontSplitBarSx(2,ii) XpBar1(ii)            YpBar1(ii)            NaN                   NaN                   0 codMatAirRot indexEle
                    XpBar1(ii)            YpBar1(ii)            XpontRadBarSx(ii)     YpontRadBarSx(ii)     NaN                   NaN                   0 codMatAirRot indexEle
                    ];
                indexEle = indexEle+1;
            end
            % external section
            if flagPosSplitSx(1,ii)==0 % left side of the rib in the central section
                if isnan(xC3k(ii)) % classical tangential rib
                    rotore = [rotore
                        XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) xxD2k(ii)             yyD2k(ii)             NaN                   NaN                   0 codMatAirRot indexEle
                        XcRibTraf2(ii)        YcRibTraf2(ii)        xxD2k(ii)             yyD2k(ii)             xpont(ii)             ypont(ii)             1 codMatAirRot indexEle
                        XcRibTraf1(ii)        YcRibTraf1(ii)        xpont(ii)             ypont(ii)             xxD1k(ii)             yyD1k(ii)             1 codMatAirRot indexEle
                        xxD1k(ii)             yyD1k(ii)             XpBar1(ii)            YpBar1(ii)            NaN                   NaN                   0 codMatAirRot indexEle
                        XpBar1(ii)            YpBar1(ii)            XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) NaN                   NaN                   0 codMatAirRot indexEle
                        xS01k(ii)             yS01k(ii)             XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    1 codMatAirRot indexEle
                        XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    NaN                   NaN                   0 codMatAirRot indexEle
                        xS02k(ii)             yS02k(ii)             XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) 1 codMatAirRot indexEle
                        ];
                    indexEle = indexEle+1;
                else % straight tangential rib
                    rotore = [rotore
                        XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) xxD2k(ii)             yyD2k(ii)             NaN                   NaN                   0 codMatAirRot indexEle
                        xC02k(ii)             yC02k(ii)             xxD2k(ii)             yyD2k(ii)             xC3k(ii)              yC3k(ii)              1 codMatAirRot indexEle
                        0                     0                     xC3k(ii)              yC3k(ii)              xC2k(ii)              yC2k(ii)              1 codMatAirRot indexEle
                        xC01k(ii)             yC01k(ii)             xC2k(ii)              yC2k(ii)              xxD1k(ii)             yyD1k(ii)             1 codMatAirRot indexEle
                        xxD1k(ii)             yyD1k(ii)             XpBar1(ii)            YpBar1(ii)            NaN                   NaN                   0 codMatAirRot indexEle
                        XpBar1(ii)            YpBar1(ii)            XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) NaN                   NaN                   0 codMatAirRot indexEle
                        xS01k(ii)             yS01k(ii)             XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    1 codMatAirRot indexEle
                        XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    NaN                   NaN                   0 codMatAirRot indexEle
                        xS02k(ii)             yS02k(ii)             XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) 1 codMatAirRot indexEle
                        ];
                    indexEle = indexEle+1;
                end
            else
                if isnan(xC3k(ii)) % classical tangential rib
                    rotore = [rotore
                        XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) xxD2k(ii)             yyD2k(ii)             NaN                   NaN                   0 codMatAirRot indexEle
                        XcRibTraf2(ii)        YcRibTraf2(ii)        xxD2k(ii)             yyD2k(ii)             xpont(ii)             ypont(ii)             1 codMatAirRot indexEle
                        XcRibTraf1(ii)        YcRibTraf1(ii)        xpont(ii)             ypont(ii)             xxD1k(ii)             yyD1k(ii)             1 codMatAirRot indexEle
                        xxD1k(ii)             yyD1k(ii)             XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) NaN                   NaN                   0 codMatAirRot indexEle
                        xS01k(ii)             yS01k(ii)             XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    1 codMatAirRot indexEle
                        XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    NaN                   NaN                   0 codMatAirRot indexEle
                        xS02k(ii)             yS02k(ii)             XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) 1 codMatAirRot indexEle
                        ];
                    indexEle = indexEle+1;
                else % straight tangential rib
                    rotore = [rotore
                        XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) xxD2k(ii)             yyD2k(ii)             NaN                   NaN                   0 codMatAirRot indexEle
                        xC02k(ii)             yC02k(ii)             xxD2k(ii)             yyD2k(ii)             xC3k(ii)              yC3k(ii)              1 codMatAirRot indexEle
                        0                     0                     xC3k(ii)              yC3k(ii)              xC2k(ii)              yC2k(ii)              1 codMatAirRot indexEle
                        xC01k(ii)             yC01k(ii)             xC2k(ii)              yC2k(ii)              xxD1k(ii)             yyD1k(ii)             1 codMatAirRot indexEle
                        xxD1k(ii)             yyD1k(ii)             XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) NaN                   NaN                   0 codMatAirRot indexEle
                        xS01k(ii)             yS01k(ii)             XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    1 codMatAirRot indexEle
                        XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    NaN                   NaN                   0 codMatAirRot indexEle
                        xS02k(ii)             yS02k(ii)             XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) 1 codMatAirRot indexEle
                        ];
                    indexEle = indexEle+1;
                end
            end
        end
    end
    % PM
    if PMdimC(ii)>0
%         Mag = [Mag
%             xPMC1b(ii) yPMC1b(ii) xPMC2b(ii) yPMC2b(ii) NaN NaN 0 codMatBar indexEle
%             xPMC2b(ii) yPMC2b(ii) xPMC2t(ii) yPMC2t(ii) NaN NaN 0 codMatBar indexEle
%             xPMC2t(ii) yPMC2t(ii) xPMC1t(ii) yPMC1t(ii) NaN NaN 0 codMatBar indexEle
%             xPMC1t(ii) yPMC1t(ii) xPMC1b(ii) yPMC1b(ii) NaN NaN 0 codMatBar indexEle
%             ];
%         indexEle = indexEle+1;
        for jj=1:PMNcC(ii)
          Mag = [Mag
            xPMC1s(jj,ii)   yPMC1s(jj,ii)    xPMC2s(jj,ii)    yPMC2s(jj,ii)   NaN NaN 0 codMatBar indexEle
            xPMC2s(jj,ii)   yPMC2s(jj,ii)    xPMC2s(jj+1,ii)  yPMC2s(jj+1,ii) NaN NaN 0 codMatBar indexEle
            xPMC2s(jj+1,ii) yPMC2s(jj+1,ii)  xPMC1s(jj+1,ii)  yPMC1s(jj+1,ii) NaN NaN 0 codMatBar indexEle
            xPMC1s(jj+1,ii) yPMC1s(jj+1,ii)  xPMC1s(jj,ii)    yPMC1s(jj,ii)   NaN NaN 0 codMatBar indexEle
            ];
        indexEle = indexEle+1;
       end
    end

    if PMdimE(ii)>0
%         Mag = [Mag
%             xPME1b(ii) yPME1b(ii) xPME2b(ii) yPME2b(ii) NaN NaN 0 codMatBar indexEle
%             xPME2b(ii) yPME2b(ii) xPME2t(ii) yPME2t(ii) NaN NaN 0 codMatBar indexEle
%             xPME2t(ii) yPME2t(ii) xPME1t(ii) yPME1t(ii) NaN NaN 0 codMatBar indexEle
%             xPME1t(ii) yPME1t(ii) xPME1b(ii) yPME1b(ii) NaN NaN 0 codMatBar indexEle
%             ];
       for jj=1:PMNcE(ii)
          Mag = [Mag
            xPME1s(jj,ii)   yPME1s(jj,ii)    xPME2s(jj,ii)    yPME2s(jj,ii)   NaN NaN 0 codMatBar indexEle
            xPME2s(jj,ii)   yPME2s(jj,ii)    xPME2s(jj+1,ii)  yPME2s(jj+1,ii) NaN NaN 0 codMatBar indexEle
            xPME2s(jj+1,ii) yPME2s(jj+1,ii)  xPME1s(jj+1,ii)  yPME1s(jj+1,ii) NaN NaN 0 codMatBar indexEle
            xPME1s(jj+1,ii) yPME1s(jj+1,ii)  xPME1s(jj,ii)    yPME1s(jj,ii)   NaN NaN 0 codMatBar indexEle
            ];
        indexEle = indexEle+1;
       end
      
    end
end

% Inserimento matrice Mag in rotore
rotore=[rotore;Mag];



