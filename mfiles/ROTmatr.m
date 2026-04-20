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

function [rotor,BLKLABELSrot,geo,mat] = ROTmatr(geo,fem,mat)

% Rotor construction.
% rotor:                  	one row per FEMM line or arc
% BLKLABELSrot.xy:          center points of FEMM blocks
% BLKLABELSrot.boundary:    one row per FEMM bounday condition
% BLKLABELSrot.BarName:     names of flux barrier blocks

p       = geo.p;
ps      = geo.ps;
th_FBS  = geo.th_FBS;
r       = geo.r;
% Ar      = geo.Ar;
lm      = geo.hc_pu*geo.g;
RotType = geo.RotType;
matFBS  = mat;
hs      = geo.hs;
g0      = geo.g;

if ~strcmp(geo.RotType,'SPM') || ~strcmp(geo.RotType,'Spoke-type') || ~strcmp(geo.RotType,'SPM-Halbach')
    mat.LayerMag.Br = mat.LayerMag.Br.*ones(1,geo.nlay);   % replicate Br in case it is scalar
end

% draw a single, straight pole or a rotor slot
if strcmp(geo.RotType,'IM')
    [geo,~,mat] = drawBar(geo,mat,fem);
else
    [geo,temp,mat] = drawPole(geo,mat,fem);
end

Ar = geo.Ar;

% initialize the matrix for geometry and labels for the total rotor (ps poles)
rotor0 = geo.rotor;
xy0    = geo.BLKLABELS.rotore.xy;

[nRow_mat,nCol_mat] = size(rotor0);
[nRow_xy,nCol_xy]   = size(xy0);

rotor     = zeros(ps*nRow_mat,nCol_mat);
BarCenter = zeros(ps*nRow_xy,nCol_xy);

% geo=geo0;

if th_FBS==0
    if strcmp(geo.RotType,'IM')
        % replicate the bar (geometry and labels)
        alphaBar = 2*pi/geo.IM.Nbars;
        for ii=1:floor((geo.IM.Nbars/(2*geo.p)*geo.ps))
            % geometry
            rotorTmp = rotateMatrix(rotor0,alphaBar*(ii-1));
            indexEle = rotorTmp(:,9);
            rotorTmp(:,9) = indexEle+max(indexEle)*(ii-1);
            rotor(1+nRow_mat*(ii-1):nRow_mat*ii,:) = rotorTmp;

            % labels
            [xtemp,ytemp] = rot_point(xy0(:,1),xy0(:,2),alphaBar*(ii-1));
            BarCenter(1+nRow_xy*(ii-1):nRow_xy*ii,:) = [xtemp,ytemp,xy0(:,3:5),0,0,xy0(:,8)];
        end
    else
        % replicate the pole (ps-1) times (geometry and labels)
        for ii=1:ps
            % geometry
            rotorTmp = rotateMatrix(rotor0,(ii-1)*pi/p);
            indexEle = rotorTmp(:,9);
            rotorTmp(:,9) = indexEle+max(indexEle)*(ii-1);
            rotor(1+nRow_mat*(ii-1):nRow_mat*ii,:) = rotorTmp;

            % labels
            if strcmp(RotType,'EESM')
                [xtemp,ytemp] = rot_point(xy0(:,1),xy0(:,2),(ii-1)*pi/p);
                magdir = atan2(xy0(:,7),xy0(:,6))+((ii-1)*pi/p-eps)+(cos((ii-2)*pi)+1)/2*pi;
                if ~mod(ii,2)
                    BarCenter(1+nRow_xy*(ii-1):nRow_xy*ii,:) = [xtemp, ytemp, xy0(:,3:5), cos(magdir), sin(magdir), -xy0(:,8)];
                else
                    BarCenter(1+nRow_xy*(ii-1):nRow_xy*ii,:) = [xtemp, ytemp, xy0(:,3:5), cos(magdir), sin(magdir), xy0(:,8)];
                end
                if ii<ps
                    BarCenter(4*ii-1,:) = NaN*ones(1,nCol_xy);
                else
                    if ps == 2*p
                        BarCenter(end-1,:) = NaN*ones(1,nCol_xy);
                    end
                end
            else
                [xtemp,ytemp] = rot_point(xy0(:,1),xy0(:,2),(ii-1)*pi/p);
                magdir = atan2(xy0(:,7),xy0(:,6))+((ii-1)*pi/p-eps)+(cos((ii-2)*pi)+1)/2*pi;
                BarCenter(1+nRow_xy*(ii-1):nRow_xy*ii,:) = [xtemp, ytemp, xy0(:,3:5), cos(magdir), sin(magdir), xy0(:,8)];
            end
        end
    end
else
    % compute ps times the shifted poles and compose the total rotor
    delta_FBS = geo.th_FBS*[-1 1];
    delta_FBS = repmat(delta_FBS,[1,ps/2]);
    geoFBS = geo;
    geoFBS.ps = 1;
    %matFBS=mat;
    for ii=1:length(delta_FBS)
        % Design of one deformed pole
        % geoFBS.delta_FBS = delta_FBS(ii);
        [rotorFBS,BLKLABELSrotFBS] = drawPoleFBS(geoFBS,matFBS,delta_FBS(ii));
        % pole rotation angle
        th_rot=pi/p*(ii-1)+sum(delta_FBS(1:ii-1))+delta_FBS(ii)/2;
        % geometry rotation
        rotorFBS = rotateMatrix(rotorFBS,th_rot);
        rotorFBS(:,9) = rotorFBS(:,9)+max(rotorFBS(:,9))*(ii-1);
        rotor(1+nRow_mat*(ii-1):nRow_mat*ii,:) = rotorFBS;
        % labels rotation
        xyFBS = BLKLABELSrotFBS.xy;
        [xtemp,ytemp] = rot_point(xyFBS(:,1),xyFBS(:,2),th_rot);
        magdir = atan2(xyFBS(:,7),xyFBS(:,6))+th_rot+pi*floor(rem(ii+1,2));
        BarCenter(1+nRow_xy*(ii-1):nRow_xy*ii,:) = [xtemp, ytemp, xy0(:,3:5), cos(magdir), sin(magdir), xy0(:,8)];
    end
end

% check the geometry matrix
[rotor] = checkPlotMatrix(rotor,1e-9);

% complete the matrix geometry (outer rotor, shaft and pole sides if needed)
if strcmp(RotType,'SPM') || strcmp(geo.RotType,'SPM-Halbach')
    re = r-lm-hs;
elseif strcmp(RotType,'EESM')
    re = geo.r-hs;
else
    re = r-hs;
end

materialCodes;
% codMatFeRot    = 5;
% codMatShaft = 7;
indexEle = max(rotor(:,9));

if (ps<2*p)
    % partial machine
    xri2 = Ar;
    yri2 = 0;
    xri3 = Ar*cos(pi/p*ps);
    yri3 = Ar*sin(pi/p*ps);
    xra2 = Ar;
    yra2 = 0;
    xra3 = Ar*cos(pi/p*ps);
    yra3 = Ar*sin(pi/p*ps);
    xre2 = re;
    yre2 = 0;
    xre3 = re*cos(pi/p*ps);
    yre3 = re*sin(pi/p*ps);
    
    if strcmp(RotType,'EESM')
        [xp3L,yp3L] = rot_point(temp.xp3,temp.yp3,-pi/p/2);
        [xp3U,yp3U] = rot_point(xp3L,yp3L,pi/p*ps);
        [xp7U,yp7U] = rot_point(temp.xp7(1),temp.yp7(1),pi/p*ps/2);
        ii = 1;
        if xre3 > xp7U
            flag = 0;
            while flag == 0
            [xp7U,yp7U] = rot_point(temp.xp7(ii),temp.yp7(ii),pi/p*ps/2);
                if (yre3-yp7U)/(xre3-xp7U) - yp7U/xp7U < tand(18) % minimum mesh angle on FEMM is 15° -> 18° for extra margin
                    ii = ii + 1;
                else
                    flag = 1;
                end
            end
        end
        [xp7L,yp7L] = rot_point(temp.xp7(ii),temp.yp7(ii),-pi/p/2);
        [xp7U,yp7U] = rot_point(xp7L,yp7L,pi/p*ps);
        % [xc4U,yc4U] = rot_point(temp.xc4,temp.yc4,pi/p*ps/2);
        % [xc4L,yc4L] = rot_point(temp.xc4,temp.yc4,-pi/p*ps/2);
        % [xc1U,yc1U] = rot_point(temp.xc1,temp.yc1,pi/p*ps/2);
        % [xc1L,yc1L] = rot_point(temp.xc1,temp.yc1,-pi/p*ps/2);
        rotor = [rotor
            0     0     xra2  yra2  NaN   NaN   0  codMatShaft   indexEle+1
            0     0     xra2  yra2  xra3  yra3 +1  codMatShaft   indexEle+1
            xra3  yra3  0     0     NaN   NaN   0  codMatShaft   indexEle+1
            xri2  yri2  xp3L  yp3L  NaN   NaN   0  codMatFeRot   indexEle+2
            xp3L -yp3L  xre2  yre2  NaN   NaN   0  codMatAirRot  indexEle+2
            xre2  yre2  xp7L -yp7L  NaN   NaN   0  codMatAirRot  indexEle+2
            xp7U  yp7U  xre3  yre3  NaN   NaN   0  codMatAirRot  indexEle+3
            xre3  yre3  xp3U  yp3U  NaN   NaN   0  codMatAirRot  indexEle+3
            xp3U  yp3U  xri3  yri3  NaN   NaN   0  codMatFeRot   indexEle+3
            ];
        if ps > 1
            for ii = 1:1:(ps-1)
            [xA,yA] = rot_point(temp.xp7(1),temp.yp7(1),pi/p/2+(ii-1)*pi/p);
            [xB,yB] = rot_point(temp.xp7(1),temp.yp7(1),pi/p/2+(ii-1)*pi/p+(1-geo.dalpha_pu)*pi/p);
            rotor = [rotor
                xA yA xB yB NaN NaN 0 codMatAirRot indexEle+3+ii
                ];
            end
        end
        if hs>0
            xre2 = r-hs;
            yre2 = 0;
            xre3 = (r-hs)*cos(pi/p*ps);
            yre3 = (r-hs)*sin(pi/p*ps);
    
            xrs2 = r;
            yrs2 = 0;
            xrs3 = (r)*cos(pi/p*ps);
            yrs3 = (r)*sin(pi/p*ps);
    
            rotor = [rotor
                xre2 yre2 xrs2 yrs2  NaN  NaN  0 codMatSleeve indexEle+4
                0    0 xrs2 yrs2 xrs3 yrs3 +1 codMatSleeve indexEle+4
                xrs3 yrs3 xre3 yre3  NaN  NaN  0 codMatSleeve indexEle+4
                0    0 xrs3 yrs3 xrs2 yrs2 -1 codMatSleeve indexEle+4
                ];
        end
    else
        rotor = [rotor
            xra2 yra2 xre2 yre2 NaN  NaN   0 codMatFeRot indexEle+1
            0    0    xre2 yre2 xre3 yre3 +1 codMatFeRot indexEle+1
            xre3 yre3 xra3 yra3 NaN  NaN   0 codMatFeRot indexEle+1
            0    0    xra3 yra3 xra2 yra2 -1 codMatFeRot indexEle+1
            0    0    xra2 yra2 NaN  NaN   0 codMatShaft indexEle+2
            0    0    xra2 yra2 xra3 yra3 +1 codMatShaft indexEle+2
            xra3 yra3 0    0    NaN  NaN   0 codMatShaft indexEle+2
            ];
        if hs>0
            xre2 = r-hs;
            yre2 = 0;
            xre3 = (r-hs)*cos(pi/p*ps);
            yre3 = (r-hs)*sin(pi/p*ps);
    
            xrs2 = r;
            yrs2 = 0;
            xrs3 = (r)*cos(pi/p*ps);
            yrs3 = (r)*sin(pi/p*ps);
    
            rotor = [rotor
                xre2 yre2 xrs2 yrs2  NaN  NaN  0 codMatSleeve indexEle+3
                0    0 xrs2 yrs2 xrs3 yrs3 +1 codMatSleeve indexEle+3
                xrs3 yrs3 xre3 yre3  NaN  NaN  0 codMatSleeve indexEle+3
                0    0 xrs3 yrs3 xrs2 yrs2 -1 codMatSleeve indexEle+3
                ];
        end
    end

    
else
    % full machine
    xre2 = re;
    yre2 = 0;
    xre3 = -re;
    yre3 = 0;

    xra2 = Ar;
    yra2 = 0;
    xra3 = -Ar;
    yra3 = 0;
    if strcmp(RotType,'EESM')
        rotor = [rotor
            0 0 xra2 yra2 xra3 yra3 1 codMatShaft indexEle+1
            0 0 xra3 yra3 xra2 yra2 1 codMatShaft indexEle+1
            ];
        for ii = 1:1:2*p
            [xA,yA] = rot_point(temp.xp7(1),temp.yp7(1),pi/p/2+ii*pi/p);
            [xB,yB] = rot_point(temp.xp7(1),temp.yp7(1),pi/p/2+ii*pi/p+(1-geo.dalpha_pu)*pi/p);
            rotor = [rotor
                xA yA xB yB NaN NaN 0 codMatAirRot indexEle+1+ii
                ];
        end
    else
        rotor = [rotor
            0 0 xre2 yre2 xre3 yre3 1 codMatFeRot indexEle+1
            0 0 xre3 yre3 xre2 yre2 1 codMatFeRot indexEle+1
            0 0 xra2 yra2 xra3 yra3 1 codMatShaft indexEle+2
            0 0 xra3 yra3 xra2 yra2 1 codMatShaft indexEle+2
            ];
    end

    if hs>0
        xrs2 = re+hs;
        yrs2 = 0;
        xrs3 = -(re+hs);
        yrs3 = 0;

        rotor = [rotor
            0 0 xrs2 yrs2 xrs3 yrs3 1 codMatSleeve indexEle+3
            0 0 xrs3 yrs3 xrs2 yrs2 1 codMatSleeve indexEle+3
            ];
    end
end

% add label for rotor iron

if strcmp(RotType,'SPM') || strcmp(geo.RotType,'SPM-Halbach')
    [xtemp,ytemp] = rot_point(mean([re Ar]),0,pi/2/p);
    BarCenter = [BarCenter; xtemp ytemp codMatFeRot,fem.res,1,NaN,NaN,NaN];
elseif strcmp(RotType,'IM')
    [xtemp,ytemp] = rot_point(mean([r-geo.IM.lt Ar]),0,pi/2/p);
    BarCenter = [BarCenter; xtemp ytemp codMatFeRot,fem.res,1,NaN,NaN,NaN];
elseif strcmp(RotType,'Spoke-type')
    [xtemp,ytemp] = rot_point(mean([r-geo.pontT-geo.PMdim(1,1) Ar]),0,pi/2/p);
    BarCenter = [BarCenter; xtemp ytemp codMatFeRot,fem.res,1,NaN,NaN,NaN];
elseif strcmp(RotType,'EESM')
    [xtemp,ytemp] = rot_point(mean([Ar Ar+geo.lyr]),0,pi/2/p);
    BarCenter = [BarCenter; xtemp ytemp codMatFeRot,fem.res,1,NaN,NaN,NaN];
else
    pointVect = sort([geo.Ar geo.B1k geo.B2k geo.r]);
    pointVect = pointVect(1:end-1)+diff(pointVect)/2;
    pointVect = pointVect(1:2:2*geo.nlay+1);
    pontFilt  = [0 fliplr(geo.pontR+geo.pontT)];
    pontFilt  = pontFilt./pontFilt;
    pointVect(1) = geo.Ar+geo.pont0/2;
    pointVect = pointVect(isnan(pontFilt));

    [xtemp,ytemp] = rot_point(pointVect',zeros(size(pointVect))',pi/2/p);
    BarCenter = [BarCenter;
        xtemp,ytemp,codMatFeRot*ones(length(xtemp),1),fem.res*ones(length(xtemp),1),ones(length(xtemp),1) nan(length(xtemp),3)];
    for ii=1:ps-1
        if length(pointVect)>1
            pointTmp = pointVect(2:end);
            [xtemp,ytemp] = rot_point(pointTmp',zeros(size(pointTmp))',pi/2/p+pi/p*(ii));
            BarCenter = [BarCenter;
                xtemp,ytemp,codMatFeRot*ones(length(xtemp),1),fem.res*ones(length(xtemp),1),ones(length(xtemp),1) nan(length(xtemp),3)];
        end
    end
end


% add label for shaft
[xtemp,ytemp] = rot_point(Ar/2,0,pi/2/p);
BarCenter = [BarCenter; xtemp ytemp codMatShaft,fem.res,1,NaN,NaN,NaN];

% add label for sleeve (if present)
if hs>0
    [xtemp,ytemp] = rot_point(r-hs/2,0,pi/2/p);
    BarCenter = [BarCenter; xtemp ytemp codMatSleeve,fem.res_traf,1,NaN,NaN,NaN];
end

% Assign label names

BarName = defineBlockNames(BarCenter,geo);

% Boundary conditions
if (ps<2*geo.p)
    codBound_periodic = 10;           % 10 = Odd or Even Periodicity
else
    codBound_periodic = -10;          % -10 = no periodicity, simulate full machine
end

% shaft boundary
[xShaftBound1,yShaftBound1] = rot_point(mean([0,Ar]),0,-90/p*pi/180);
[xShaftBound2,yShaftBound2] = rot_point(mean([0,Ar]),0,(ps-1/2)*180/p*pi/180);
% rotor and sleeve boundary
if strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'SPM-Halbach')
    [xRotBound1,yRotBound1] = rot_point(mean([Ar,r-hs]),0,-90/p*pi/180);
    [xRotBound2,yRotBound2] = rot_point(mean([Ar,r-hs]),0,(ps-1/2)*180/p*pi/180);
    xRotAirBound1 = NaN; 
    yRotAirBound1 = NaN;
    xRotAirBound2 = NaN; 
    yRotAirBound2 = NaN;
    [xSleeveBound1,ySleeveBound1] = rot_point(mean([r-hs,r]),0,-90/p*pi/180);
    [xSleeveBound2,ySleeveBound2] = rot_point(mean([r-hs,r]),0,(ps-1/2)*180/p*pi/180);
elseif strcmp(geo.RotType,'EESM')
    [xRotBound1,yRotBound1] = rot_point(mean([Ar,Ar+geo.lyr]),0,-90/p*pi/180);
    [xRotBound2,yRotBound2] = rot_point(mean([Ar,Ar+geo.lyr]),0,(ps-1/2)*180/p*pi/180);
    [xRotAirBound1,yRotAirBound1] = rot_point(mean([Ar+geo.lyr,r-hs]),0,-90/p*pi/180);
    [xRotAirBound2,yRotAirBound2] = rot_point(mean([Ar+geo.lyr,r-hs]),0,(ps-1/2)*180/p*pi/180);
    [xSleeveBound1,ySleeveBound1] = rot_point(mean([r-hs,r]),0,-90/p*pi/180);
    [xSleeveBound2,ySleeveBound2] = rot_point(mean([r-hs,r]),0,(ps-1/2)*180/p*pi/180);
else
    [xRotBound1,yRotBound1] = rot_point(mean([Ar,r-lm-hs]),0,-90/p*pi/180);          % for SPM motor
    [xRotBound2,yRotBound2] = rot_point(mean([Ar,r-lm-hs]),0,(ps-1/2)*180/p*pi/180);
    xRotAirBound1 = NaN; 
    yRotAirBound1 = NaN;
    xRotAirBound2 = NaN; 
    yRotAirBound2 = NaN;
    [xSleeveBound1,ySleeveBound1] = rot_point(mean([r-lm-hs,r]),0,-90/p*pi/180);
    [xSleeveBound2,ySleeveBound2] = rot_point(mean([r-lm-hs,r]),0,(ps-1/2)*180/p*pi/180);
end


%%% OUTPUT DATA %%%
%%%%%%%%%%%%%%%%%%%

geo.rotor = rotor;

%%% Block centers %%%
BLKLABELSrot.xy     =   BarCenter;
BLKLABELSrot.BarName =   BarName';

% Boundaries %%%
BLKLABELSrot.boundary = [
    xShaftBound1  yShaftBound1  codBound_periodic;
    xShaftBound2  yShaftBound2  codBound_periodic;
    xRotBound1    yRotBound1    codBound_periodic;
    xRotBound2    yRotBound2    codBound_periodic;
    xRotAirBound1 yRotAirBound1 codBound_periodic;
    xRotAirBound2 yRotAirBound2 codBound_periodic;
    xSleeveBound1 ySleeveBound1 codBound_periodic;
    xSleeveBound2 ySleeveBound2 codBound_periodic;
    ];

BLKLABELSrot.boundary = BLKLABELSrot.boundary(~isnan(BLKLABELSrot.boundary(:,1)),:);    %remove NaN rows 

% Rotate boundary selection points
[xtemp,ytemp]=rot_point(BLKLABELSrot.boundary(:,1),BLKLABELSrot.boundary(:,2),90/p*pi/180);
BLKLABELSrot.boundary=[xtemp,ytemp,BLKLABELSrot.boundary(:,3:end)];
%

