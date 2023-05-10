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
lm      = geo.lm;
RotType = geo.RotType;
matFBS  = mat;
hs      = geo.hs;

if ~strcmp(geo.RotType,'SPM')
    mat.LayerMag.Br = mat.LayerMag.Br.*ones(1,geo.nlay);   % replicate Br in case it is scalar
end

% draw a single, straight pole or a rotor slot
if strcmp(geo.RotType,'IM')
    [geo,~,mat] = drawBar(geo,mat,fem);
else
    [geo,~,mat] = drawPole(geo,mat,fem);
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
            [xtemp,ytemp] = rot_point(xy0(:,1),xy0(:,2),(ii-1)*pi/p);
            magdir = atan2(xy0(:,7),xy0(:,6))+((ii-1)*pi/p-eps)+(cos((ii-2)*pi)+1)/2*pi;
            BarCenter(1+nRow_xy*(ii-1):nRow_xy*ii,:) = [xtemp, ytemp, xy0(:,3:5), cos(magdir), sin(magdir), xy0(:,8)];
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
        geoFBS.delta_FBS = delta_FBS(ii);
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
if strcmp(RotType,'SPM')
    re = r-lm;
else
    re = r;
end

materialCodes;
% codMatFeRot    = 5;
% codMatShaft = 7;
indexEle = max(rotor(:,9));

if (ps<2*p)
    % partial machine
    xre2 = re;
    yre2 = 0;
    xre3 = re*cos(pi/p*ps);
    yre3 = re*sin(pi/p*ps);
    
    xra2 = Ar;
    yra2 = 0;
    xra3 = Ar*cos(pi/p*ps);
    yra3 = Ar*sin(pi/p*ps);

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
        xrs2 = re+hs;
        yrs2 = 0;
        xrs3 = (re+hs)*cos(pi/p*ps);
        yrs3 = (re+hs)*sin(pi/p*ps);

        rotor = [rotor
            xre2 yre2 xrs2 yrs2  NaN  NaN  0 codMatSleeve indexEle+3
               0    0 xrs2 yrs2 xrs3 yrs3 +1 codMatSleeve indexEle+3
            xrs3 yrs3 xre3 yre3  NaN  NaN  0 codMatSleeve indexEle+3
               0    0 xrs3 yrs3 xrs2 yrs2 -1 codMatSleeve indexEle+3
            ];
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
    
    rotor = [rotor
        0 0 xre2 yre2 xre3 yre3 1 codMatFeRot indexEle+1
        0 0 xre3 yre3 xre2 yre2 1 codMatFeRot indexEle+1
        0 0 xra2 yra2 xra3 yra3 1 codMatShaft indexEle+2
        0 0 xra3 yra3 xra2 yra2 1 codMatShaft indexEle+2
        ];

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

if strcmp(RotType,'SPM')
    [xtemp,ytemp] = rot_point(mean([re Ar]),0,pi/2/p);
    BarCenter = [BarCenter; xtemp ytemp codMatFeRot,fem.res,1,NaN,NaN,NaN];
elseif strcmp(RotType,'IM')
    [xtemp,ytemp] = rot_point(mean([r-geo.IM.lt Ar]),0,pi/2/p);
    BarCenter = [BarCenter; xtemp ytemp codMatFeRot,fem.res,1,NaN,NaN,NaN];
else
    pointVect = sort([geo.Ar geo.B1k geo.B2k geo.r]);
    pointVect = pointVect(1:end-1)+diff(pointVect)/2;
    pointVect = pointVect(1:2:2*geo.nlay+1);
    pontFilt  = [0 fliplr(geo.pontR+geo.pontT)];
    pontFilt  = pontFilt./pontFilt;
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
    [xtemp,ytemp] = rot_point(r+hs/2,0,pi/2/p);
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
% rotor boundary
if strcmp(geo.RotType,'SPM')
    [xRotBound1,yRotBound1] = rot_point(mean([Ar,r+hs]),0,-90/p*pi/180);
    [xRotBound2,yRotBound2] = rot_point(mean([Ar,r+hs]),0,(ps-1/2)*180/p*pi/180);
else
    [xRotBound1,yRotBound1] = rot_point(mean([Ar,r-lm+hs]),0,-90/p*pi/180);          % for SPM motor
    [xRotBound2,yRotBound2] = rot_point(mean([Ar,r-lm+hs]),0,(ps-1/2)*180/p*pi/180);
end
% sleeve boundary


%%% OUTPUT DATA %%%
%%%%%%%%%%%%%%%%%%%

geo.rotor = rotor;

%%% Block centers %%%
BLKLABELSrot.xy     =   BarCenter;
BLKLABELSrot.BarName =   BarName';

% Boundaries %%%
BLKLABELSrot.boundary = [
    xShaftBound1 yShaftBound1 codBound_periodic;
    xShaftBound2 yShaftBound2 codBound_periodic;
    xRotBound1   yRotBound1   codBound_periodic;
    xRotBound2   yRotBound2   codBound_periodic
    ];
% Rotate boundary selection points
[xtemp,ytemp]=rot_point(BLKLABELSrot.boundary(:,1),BLKLABELSrot.boundary(:,2),90/p*pi/180);
BLKLABELSrot.boundary=[xtemp,ytemp,BLKLABELSrot.boundary(:,3:end)];
%

