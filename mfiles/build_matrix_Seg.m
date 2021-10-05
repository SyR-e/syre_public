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

flag_segV       = temp.flag_segV;
splitRib        = temp.splitRib;

if any(splitRib)
    flag_cent       = temp.flag_cent;
    flag_ext        = temp.flag_ext;
    flag_shift      = temp.flag_shift;
    flag_shiftUP    = temp.flag_shiftUP;
    xI01k=temp.xI01k;
    yI01k=temp.yI01k;
    xI02k=temp.xI02k;
    yI02k=temp.yI02k;

    
end
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

if strcmp(geo.RotType,'Seg')
        xC1k=temp.xC1k;
        yC1k=temp.yC1k;
        xC2k=temp.xC2k;
        yC2k=temp.yC2k;
        xC3k=temp.xC3k;
        yC3k=temp.yC3k;
        xC4k=temp.xC4k;
        yC4k=temp.yC4k;
        xC01k=temp.xC01k;
        yC01k=temp.yC01k;
        xC02k=temp.xC02k;
        yC02k=temp.yC02k;
        
        
        %%%%%%%
        %         XpontSplitBarSx(1,:)=temp.XpontSplitBarSx(1,:);
        %         YpontSplitBarSx(1,:)=temp.YpontSplitBarSx(1,:);
        %         XpontSplitDx(1,:)=temp.XpontSplitDx(1,:);
        %         YpontSplitDx(1,:)=temp.YpontSplitDx(1,:);
        %         XpontSplitBarDx(1,:)=temp.XpontSplitBarDx(1,:);
        %         YpontSplitBarDx(1,:)=temp.YpontSplitBarDx(1,:);
        %         XpontSplitDx(1,:)=temp.XpontSplitDx(1,:);
        %         YpontSplitDx(1,:)=temp.YpontSplitDx(1,:);

        %
        %         XpontSplitBarSx(2,:)=temp.XpontSplitBarSx(2,:);
        %         YpontSplitBarSx(2,:)=temp.YpontSplitBarSx(2,:);
        %         XpontSplitSx(2,:)=temp.XpontSplitSx(2,:);
        %         YpontSplitSx(2,:)=temp.YpontSplitSx(2,:);
        %         XpontSplitBarDx(2,:)=temp.XpontSplitBarDx(2,:);
        %         YpontSplitBarDx(2,:)=temp.YpontSplitBarDx(2,:);
        %         XpontSplitDx(2,:)=temp.XpontSplitDx(2,:);
        %         YpontSplitDx(2,:)=temp.YpontSplitDx(2,:);
end

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
    %% first section: from q-axis to the arms. counterclockwise contour from (XpBar1,YpBar1) to (XpBar2,YpBar2)
    if pontR(ii)~=0
        if splitRib(ii)
%             if flag_pont(ii)==1
                if (flag_cent(ii)==1)           % end-rib goes on the central barrier 
                    rotore = [rotore
                        XpontRadBarSx(ii)     YpontRadBarSx(ii)     XpontRadBarDx(ii)     YpontRadBarDx(ii)     NaN                     NaN                   0 codMatAirRot indexEle
                        XpontRadBarDx(ii)     YpontRadBarDx(ii)     XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) NaN                     NaN                   0 codMatAirRot indexEle
                        xI02k(ii)             yI02k(ii)             XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) XpontSplitDx(2,ii)      YpontSplitDx(2,ii)    1 codMatAirRot indexEle
                        XpontSplitDx(2,ii)    YpontSplitDx(2,ii)    XpontSplitSx(2,ii)    YpontSplitSx(2,ii)    NaN                     NaN                   0 codMatAirRot indexEle
                        xI01k(ii)             yI01k(ii)             XpontSplitSx(2,ii)    YpontSplitSx(2,ii)    XpontSplitBarSx(2,ii)   YpontSplitBarSx(2,ii) 1 codMatAirRot indexEle
                        XpontSplitBarSx(2,ii) YpontSplitBarSx(2,ii) XpontRadBarSx(ii)     YpontRadBarSx(ii)     NaN                     NaN                   0 codMatAirRot indexEle
                        ];
                    indexEle = indexEle+1;
                end
                %                  if (flag_cent(ii)==1 && flag_ext(ii)==1)
                %                 rotore = [rotore
                %                     XpontRadBarSx(ii)     YpontRadBarSx(ii)     XpontRadBarDx(ii)     YpontRadBarDx(ii)     NaN NaN 0 codMatAirRot indexEle
                %                     XpontRadBarDx(ii)     YpontRadBarDx(ii)     XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) NaN NaN 0 codMatAirRot indexEle
                %                     xI02k(ii)             yI02k(ii)             XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) XpontSplitDx(2,ii)   YpontSplitDx(2,ii) 1 codMatAirRot indexEle
                %                     XpontSplitSx(2,ii)    YpontSplitSx(2,ii)    XpontSplitDx(2,ii)    YpontSplitDx(2,ii)    NaN NaN 0 codMatAirRot indexEle
                %                     xI01k(ii)             yI01k(ii)             XpontSplitSx(2,ii)    YpontSplitSx(2,ii)    XpontSplitBarSx(2,ii)   YpontSplitBarSx(2,ii) 1 codMatAirRot indexEle
                %                     XpontSplitBarSx(2,ii) YpontSplitBarSx(2,ii) XpontRadBarSx(ii)     YpontRadBarSx(ii)     NaN NaN 0 codMatAirRot indexEle
                %                     ];
                %                 indexEle = indexEle+1;
                %                  end
                if (flag_cent(ii)==0 && flag_ext(ii)==1)              % end-rib goes on the external barrier
                    if flag_shiftUP(ii)==1                            % rib with positive offset (shifted towards the exteranal barrier)
                        rotore = [rotore
                            XpontRadBarSx(ii)     YpontRadBarSx(ii)     XpontRadBarDx(ii)     YpontRadBarDx(ii)     NaN NaN 0 codMatAirRot indexEle
                            XpontRadBarDx(ii)     YpontRadBarDx(ii)     XpBar2(ii)            YpBar2(ii)            NaN NaN 0 codMatAirRot indexEle
                            XpBar2(ii)            YpBar2(ii)            XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) NaN NaN 0 codMatAirRot indexEle
                            ];
                    else
                        rotore = [rotore
                            XpontRadBarSx(ii)     YpontRadBarSx(ii)     XpontRadBarDx(ii)     YpontRadBarDx(ii)     NaN NaN 0 codMatAirRot indexEle
                            XpontRadBarDx(ii)     YpontRadBarDx(ii)     XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) NaN NaN 0 codMatAirRot indexEle
                            ];
                    end
                    rotore = [rotore
                        xI02k(ii)             yI02k(ii)             XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) XpontSplitDx(2,ii)    YpontSplitDx(2,ii)    1 codMatAirRot indexEle
                        XpontSplitDx(2,ii)    YpontSplitDx(2,ii)    XpontSplitSx(2,ii)    YpontSplitSx(2,ii)    NaN                   NaN                   0 codMatAirRot indexEle
                        xI01k(ii)             yI01k(ii)             XpontSplitSx(2,ii)    YpontSplitSx(2,ii)    XpontSplitBarSx(2,ii) YpontSplitBarSx(2,ii) 1 codMatAirRot indexEle
                        XpontSplitBarSx(2,ii) YpontSplitBarSx(2,ii) XpBar1(ii)            YpBar1(ii)            NaN                   NaN                   0 codMatAirRot indexEle
                        XpBar1(ii)            YpBar1(ii)            XpontRadBarSx(ii)     YpontRadBarSx(ii)     NaN                   NaN                   0 codMatAirRot indexEle
                        ];
                    indexEle = indexEle+1;
                end
                
                % end
                if ((isfinite(geo.RotorFilletTan1(ii))) && (strcmp(geo.RotType,'Seg')))

                    if flag_ext(ii)==1
                        rotore = [rotore
                            xC1k(ii)              yC1k(ii)              XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) NaN                   NaN                   0 codMatAirRot indexEle
                            xS01k(ii)             yS01k(ii)             XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    1 codMatAirRot indexEle
                            XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    NaN                   NaN                   0 codMatAirRot indexEle
                            xS02k(ii)             yS02k(ii)             XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) 1 codMatAirRot indexEle
                            ];
                        
                    else
                        rotore = [rotore
                            xC1k(ii)              yC1k(ii)              XpBar1(ii)            YpBar1(ii)            NaN                   NaN                   0 codMatAirRot indexEle
                            XpBar1(ii)            YpBar1(ii)            XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) NaN                   NaN                   0 codMatAirRot indexEle
                            xS01k(ii)             yS01k(ii)             XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    1 codMatAirRot indexEle
                            XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    NaN                   NaN                   0 codMatAirRot indexEle
                            xS02k(ii)             yS02k(ii)             XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) 1 codMatAirRot indexEle
                            ];
                    end
                    
                    if flag_shift(ii)==1                                        % rib with negative offset (shifted towards the center)
                        rotore = [rotore
                            XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) XpBar2(ii)  YpBar2(ii)            NaN                   NaN                   0 codMatAirRot indexEle
                            XpBar2(ii)            YpBar2(ii)            xC4k(ii)    yC4k(ii)              NaN                   NaN                   0 codMatAirRot indexEle
                            ];
                    else
                        rotore = [rotore
                            XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) xC4k(ii)              yC4k(ii)              NaN                   NaN                   0 codMatAirRot indexEle
                            ];
                        
                    end
                else
 
                    if flag_ext(ii)==1                              % end-rib goes on the external barrier
                        rotore = [rotore
                            xxD1k(ii)             yyD1k(ii)             XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) NaN                   NaN                   0 codMatAirRot indexEle
                            xS01k(ii)             yS01k(ii)             XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    1 codMatAirRot indexEle
                            XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    NaN                   NaN                   0 codMatAirRot indexEle
                            xS02k(ii)             yS02k(ii)             XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) 1 codMatAirRot indexEle
                            ];
                    else
                        rotore = [rotore
                            xxD1k(ii)             yyD1k(ii)             XpBar1(ii)            YpBar1(ii)            NaN                   NaN                   0 codMatAirRot indexEle
                            XpBar1(ii)            YpBar1(ii)            XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) NaN                   NaN                   0 codMatAirRot indexEle
                            xS01k(ii)             yS01k(ii)             XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    1 codMatAirRot indexEle
                            XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    NaN                   NaN                   0 codMatAirRot indexEle
                            xS02k(ii)             yS02k(ii)             XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) 1 codMatAirRot indexEle
                            ];
                    end
                    
                    if flag_shift(ii)==1                                            % rib with negative offset (shifted towards the center)
                        rotore = [rotore
                            XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) XpBar2(ii)  YpBar2(ii)            NaN                   NaN                   0 codMatAirRot indexEle
                            XpBar2(ii)            YpBar2(ii)            xxD2k(ii)   yyD2k(ii)              NaN                   NaN                   0 codMatAirRot indexEle
                            ];
                    else
                        rotore = [rotore
                            XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) xxD2k(ii)              yyD2k(ii)              NaN                   NaN                   0 codMatAirRot indexEle
                            ];
                    end
                    

                end
%             else
%                 rotore = [rotore
%                     XpontRadBarSx(ii)     YpontRadBarSx(ii)     XpontRadBarDx(ii)     YpontRadBarDx(ii)     NaN NaN 0 codMatAirRot indexEle
%                     XpontRadBarDx(ii)     YpontRadBarDx(ii)     XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) NaN NaN 0 codMatAirRot indexEle
%                     XpontSplitBarDx(2,ii) YpontSplitBarDx(2,ii) XpontSplitDx(2,ii)    YpontSplitDx(2,ii)    NaN NaN 0 codMatAirRot indexEle
%                     XpontSplitDx(2,ii)    YpontSplitDx(2,ii)    XpontSplitSx(2,ii)    YpontSplitSx(2,ii)    NaN NaN 0 codMatAirRot indexEle
%                     XpontSplitSx(2,ii)    YpontSplitSx(2,ii)    XpontSplitBarSx(2,ii) YpontSplitBarSx(2,ii) NaN NaN 0 codMatAirRot indexEle
%                     XpontSplitBarSx(2,ii) YpontSplitBarSx(2,ii) XpontRadBarSx(ii)     YpontRadBarSx(ii)     NaN NaN 0 codMatAirRot indexEle
%                     ];
%                 indexEle = indexEle+1;
%                 
%                 rotore = [rotore
%                     XpBar1(ii)            YpBar1(ii)            XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) NaN NaN 0 codMatAirRot indexEle
%                     XpontSplitBarSx(1,ii) YpontSplitBarSx(1,ii) XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    NaN NaN 0 codMatAirRot indexEle
%                     XpontSplitSx(1,ii)    YpontSplitSx(1,ii)    XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    NaN NaN 0 codMatAirRot indexEle
%                     XpontSplitDx(1,ii)    YpontSplitDx(1,ii)    XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) NaN NaN 0 codMatAirRot indexEle
%                     XpontSplitBarDx(1,ii) YpontSplitBarDx(1,ii) XpBar2(ii)            YpBar2(ii)            NaN NaN 0 codMatAirRot indexEle
%                     ];
%             end
        elseif flag_segV(ii)==0
            rotore = [rotore
                XpBar1(ii)        YpBar1(ii)        XpontRadBarSx(ii) YpontRadBarSx(ii) NaN NaN 0 codMatAirRot indexEle
                xS01k(ii)          yS01k(ii)          XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadSx(ii)    YpontRadSx(ii)     1 codMatAirRot indexEle
                %XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadSx(ii)    YpontRadSx(ii)    NaN NaN 0 codMatAirRot indexEle
                XpontRadSx(ii)    YpontRadSx(ii)    XpontRadDx(ii)    YpontRadDx(ii)    NaN NaN 0 codMatAirRot indexEle
                xS02k(ii)          yS02k(ii)          XpontRadDx(ii)    YpontRadDx(ii)    XpontRadBarDx(ii) YpontRadBarDx(ii) 1 codMatAirRot indexEle
                %XpontRadDx(ii)    YpontRadDx(ii)    XpontRadBarDx(ii) YpontRadBarDx(ii) NaN NaN 0 codMatAirRot indexEle
                XpontRadBarDx(ii) YpontRadBarDx(ii) XpBar2(ii)        YpBar2(ii)        NaN NaN 0 codMatAirRot indexEle
                ];
        end
    elseif flag_segV(ii)==0
        rotore = [rotore
            XpBar1(ii)        YpBar1(ii)        XpontRadBarSx(ii) YpontRadBarSx(ii) NaN NaN 0 codMatAirRot indexEle
            XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadBarDx(ii) YpontRadBarDx(ii) NaN NaN 0 codMatAirRot indexEle
            XpontRadBarDx(ii) YpontRadBarDx(ii) XpBar2(ii)        YpBar2(ii)        NaN NaN 0 codMatAirRot indexEle
            ];
    end
    
    
    %% second section of the layer: arms and barrier edges, from (XpBar2,YpBar2) to (XpBar1,YpBar1)
    if ((isfinite(geo.RotorFilletTan1(ii))) && (strcmp(geo.RotType,'Seg')))
        if pontR(ii)~=0
             if splitRib(ii)==1
            rotore = [rotore
                %XpontSplitBarDx(1,ii)   YpontSplitBarDx(1,ii)   xC4k(ii)                yC4k(ii)                NaN      NaN      0 codMatAirRot indexEle
                xC02k(ii)               yC02k(ii)               xC4k(ii)                yC4k(ii)                xC3k(ii) yC3k(ii) 1 codMatAirRot indexEle
                xC3k(ii)                yC3k(ii)                xC2k(ii)                yC2k(ii)                NaN      NaN      0 codMatAirRot indexEle
                xC01k(ii)               yC01k(ii)               xC2k(ii)                yC2k(ii)                xC1k(ii) yC1k(ii) 1 codMatAirRot indexEle
                %xC1k(ii)                yC1k(ii)                XpontSplitBarSx(1,ii)   YpontSplitBarSx(1,ii)   NaN      NaN      0 codMatAirRot indexEle
                ];
             elseif  flag_segV(ii)==0 
            rotore = [rotore
                XpBar2(ii)              YpBar2(1,ii)            xC4k(ii)                yC4k(ii)                NaN      NaN      0 codMatAirRot indexEle
                xC02k(ii)               yC02k(ii)               xC4k(ii)                yC4k(ii)                xC3k(ii) yC3k(ii) 1 codMatAirRot indexEle
                xC3k(ii)                yC3k(ii)                xC2k(ii)                yC2k(ii)                NaN      NaN      0 codMatAirRot indexEle
                xC01k(ii)               yC01k(ii)               xC2k(ii)                yC2k(ii)                xC1k(ii) yC1k(ii) 1 codMatAirRot indexEle
                xC1k(ii)                yC1k(ii)                XpBar1(ii)              YpBar1(ii)   NaN      NaN      0 codMatAirRot indexEle
                ];
            elseif flag_segV(ii)==1
           rotore = [rotore
                %XpontRadDx(ii)    YpontRadDx(ii)    XpontRadBarDx(ii) YpontRadBarDx(ii) NaN NaN 0 codMatAirRot indexEle
                XpontRadBarDx(ii) YpontRadBarDx(ii)  xC4k(ii)                yC4k(ii)         NaN NaN 0 codMatAirRot indexEle
                xC02k(ii)               yC02k(ii)               xC4k(ii)                yC4k(ii)                xC3k(ii) yC3k(ii) 1 codMatAirRot indexEle
                xC3k(ii)                yC3k(ii)                xC2k(ii)                yC2k(ii)                NaN      NaN      0 codMatAirRot indexEle
                xC01k(ii)               yC01k(ii)               xC2k(ii)                yC2k(ii)                xC1k(ii) yC1k(ii) 1 codMatAirRot indexEle
                xC1k(ii) yC1k(ii)         XpontRadBarSx(ii) YpontRadBarSx(ii) NaN NaN 0 codMatAirRot indexEle
                xS01k(ii)          yS01k(ii)         XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadSx(ii)    YpontRadSx(ii)     1 codMatAirRot indexEle
                %XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadSx(ii)    YpontRadSx(ii)    NaN NaN 0 codMatAirRot indexEle
                XpontRadSx(ii)    YpontRadSx(ii)     XpontRadDx(ii)    YpontRadDx(ii)    NaN NaN 0 codMatAirRot indexEle
                xS02k(ii)          yS02k(ii)         XpontRadDx(ii)    YpontRadDx(ii)    XpontRadBarDx(ii) YpontRadBarDx(ii) 1 codMatAirRot indexEle
                ];
             end
        else
            rotore = [rotore
                B1k(ii)    0          B2k(ii)    0          NaN      NaN      0 codMatAirRot indexEle
                XpBar2(ii) YpBar2(ii) xC4k(ii)   yC4k(ii)   NaN      NaN      0 codMatAirRot indexEle
                xC02k(ii)  yC02k(ii)  xC4k(ii)   yC4k(ii)   xC3k(ii) yC3k(ii) 1 codMatAirRot indexEle
                xC3k(ii)   yC3k(ii)   xC2k(ii)   yC2k(ii)   NaN      NaN      0 codMatAirRot indexEle
                xC01k(ii)  yC01k(ii)  xC2k(ii)   yC2k(ii)   xC1k(ii) yC1k(ii) 1 codMatAirRot indexEle
                xC1k(ii)   yC1k(ii)   XpBar1(ii) YpBar1(ii) NaN      NaN      0 codMatAirRot indexEle
                ];
        end
    else
        if splitRib(ii)==1
            rotore = [rotore
                %XpontSplitBarDx(1,ii)       YpontSplitBarDx(1,ii)       xxD2k(ii)               yyD2k(ii)             NaN       NaN       0 codMatAirRot indexEle
                XcRibTraf2(ii)              YcRibTraf2(ii)              xxD2k(ii)               yyD2k(ii)             xpont(ii) ypont(ii) 1 codMatAirRot indexEle
                XcRibTraf1(ii)              YcRibTraf1(ii)              xpont(ii)               ypont(ii)             xxD1k(ii) yyD1k(ii) 1 codMatAirRot indexEle
                %xxD1k(ii)                   yyD1k(ii)                   XpontSplitBarSx(1,ii)   YpontSplitBarSx(1,ii) NaN       NaN       0 codMatAirRot indexEle
                ];
        elseif flag_segV(ii)==0
            rotore = [rotore
                XpBar2(ii)     YpBar2(ii)     xxD2k(ii)  yyD2k(ii)  NaN       NaN       0 codMatAirRot indexEle
                XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii)  yyD2k(ii)  xpont(ii) ypont(ii) 1 codMatAirRot indexEle
                XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii)  ypont(ii)  xxD1k(ii) yyD1k(ii) 1 codMatAirRot indexEle
                xxD1k(ii)      yyD1k(ii)      XpBar1(ii) YpBar1(ii) NaN       NaN       0 codMatAirRot indexEle
                ];
        elseif flag_segV(ii)==1
            if pontR(ii)~=0
                rotore = [rotore
                    %XpontRadDx(ii)    YpontRadDx(ii)    XpontRadBarDx(ii) YpontRadBarDx(ii) NaN NaN 0 codMatAirRot indexEle
                    XpontRadBarDx(ii) YpontRadBarDx(ii)  xxD2k(ii)         yyD2k(ii)         NaN NaN 0 codMatAirRot indexEle
                    XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii)  yyD2k(ii)  xpont(ii) ypont(ii) 1 codMatAirRot indexEle
                    XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii)  ypont(ii)  xxD1k(ii) yyD1k(ii) 1 codMatAirRot indexEle
                    xxD1k(ii)          yyD1k(ii)         XpontRadBarSx(ii) YpontRadBarSx(ii) NaN NaN 0 codMatAirRot indexEle
                    xS01k(ii)          yS01k(ii)         XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadSx(ii)    YpontRadSx(ii)     1 codMatAirRot indexEle
                    %XpontRadBarSx(ii) YpontRadBarSx(ii) XpontRadSx(ii)    YpontRadSx(ii)    NaN NaN 0 codMatAirRot indexEle
                    XpontRadSx(ii)    YpontRadSx(ii)     XpontRadDx(ii)    YpontRadDx(ii)    NaN NaN 0 codMatAirRot indexEle
                    xS02k(ii)          yS02k(ii)         XpontRadDx(ii)    YpontRadDx(ii)    XpontRadBarDx(ii) YpontRadBarDx(ii) 1 codMatAirRot indexEle
                    ];
            else
                rotore = [rotore
                    B1k(ii)    0          B2k(ii)    0          NaN      NaN      0 codMatAirRot indexEle
                    XpBar2(ii)     YpBar2(ii)     xxD2k(ii)  yyD2k(ii)  NaN       NaN       0 codMatAirRot indexEle
                    XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii)  yyD2k(ii)  xpont(ii) ypont(ii) 1 codMatAirRot indexEle
                    XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii)  ypont(ii)  xxD1k(ii) yyD1k(ii) 1 codMatAirRot indexEle
                    xxD1k(ii)      yyD1k(ii)      XpBar1(ii) YpBar1(ii) NaN       NaN       0 codMatAirRot indexEle
                    ];
            end
        end
    end
    indexEle = indexEle+1;
    
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



