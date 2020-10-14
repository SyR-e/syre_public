% Copyright 2019
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [geo,mat,temp] = PMdefinition_Seg(geo,mat,temp)
% 
% [geo,temp] = PMdefinition_Seg(geo,temp,mat)
%

% barrier arms
XpBar1 = temp.XpBar1;
YpBar1 = temp.YpBar1;
xxD1k  = temp.xxD1k;
yyD1k  = temp.yyD1k;
xxD2k  = temp.xxD2k;
yyD2k  = temp.yyD2k;
XpBar2 = temp.XpBar2;
YpBar2 = temp.YpBar2;
xpont  = temp.xpont;
ypont  = temp.ypont;
% B1k    = temp.B1k;
% B2k    = temp.B2k;

% Rad Rib
YpontRadSx    = temp.YpontRadSx;
XpontRadSx    = temp.XpontRadSx;
YpontRadDx    = temp.YpontRadDx;
XpontRadDx    = temp.XpontRadDx;
XpontRadBarDx = temp.XpontRadBarDx;
YpontRadBarDx = temp.YpontRadBarDx;
XpontRadBarSx = temp.XpontRadBarSx;
YpontRadBarSx = temp.YpontRadBarSx;

% Split Rad Rib
XpontSplitBarSx = temp.XpontSplitBarSx;
YpontSplitBarSx = temp.YpontSplitBarSx;
XpontSplitBarDx = temp.XpontSplitBarDx;
YpontSplitBarDx = temp.YpontSplitBarDx;
XpontSplitDx    = temp.XpontSplitDx;
YpontSplitDx    = temp.YpontSplitDx;
XpontSplitSx    = temp.XpontSplitSx;
YpontSplitSx    = temp.YpontSplitSx;

% from geo
p         = geo.p;
nlay      = geo.nlay;
pont0     = geo.pont0;
hc        = geo.hc;
FBS       = abs(geo.delta_FBS);
PMdimC    = geo.PMdim(1,:);
PMdimE    = geo.PMdim(2,:);
% betaC     = geo.PMdim(1,:);
% betaE     = geo.PMdim(2,:);
flagPMFBS = geo.flagPMFBS;
deltaFBS  = geo.delta_FBS;
PMclearC  = geo.PMclear(1,:);
PMclearE  = geo.PMclear(2,:);

% initialize variables for PM area.
%   C --> central portion of the barrier
%   E --> external portion of the barrier (arms)

% initialize points for PM definition
%   1 --> inner side of the barrier (1==sx)
%   2 --> outer side of the barrier (2==dx)
%   t --> top side of the magnet
%   b --> bottom side of the magnet

% point definition - central section

b = YpontSplitBarSx(2,:)-YpontRadBarSx-pont0/4;
AreaC = hc.*b;
if FBS>0
    if flagPMFBS
        AreaCMax = AreaC;
    else
        AreaCMax = geo.AreaCMax;
    end
else
    AreaCMax = AreaC;
end

if min(PMdimC)<0 % per-unit input (during optimization)
    PMdimC = -PMdimC.*b;
end

PMdimC(PMdimC>b)=b(PMdimC>b);
PMdimC(PMdimC<0)=0;
PMdimC=floor(PMdimC*100)/100;

AreaC = PMdimC.*hc;

% Check clearance size: the minimum PM thickness is set to pont0
PMclearLim = hc-pont0*ones(size(hc));
PMclearC(PMclearC>PMclearLim) = PMclearLim(PMclearC>PMclearLim);
PMclearE(PMclearE>PMclearLim) = PMclearLim(PMclearE>PMclearLim);

xPMC1b = XpontRadBarSx+PMclearC;
yPMC1b = YpontRadBarSx;
xPMC2b = XpontRadBarDx;
yPMC2b = YpontRadBarDx;

xPMC2t = xPMC2b;
yPMC2t = yPMC2b+PMdimC;
xPMC1t = xPMC1b;
yPMC1t = yPMC2t;

% point definition - external section
xPME2b = XpBar2+pont0*cos(pi/2/p+deltaFBS/2);
yPME2b = YpBar2+pont0*sin(pi/2/p+deltaFBS/2);
[a1,b1,c1] = retta_per_2pti(XpBar1,YpBar1,xxD1k,yyD1k);
a2 = tan(pi/2/p+deltaFBS/2-pi/2).*ones(1,nlay);
b2 = -1*ones(1,nlay);
c2 = YpBar2+pont0*sin(pi/2/p+deltaFBS/2)-a2.*(XpBar2+pont0*cos(pi/2/p+deltaFBS/2));
[xPME1b,yPME1b]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);

filt=YpBar1>yPME1b;
if sum(filt)
    xPME1b(filt) = XpBar1(filt)+pont0*cos(pi/2/p+deltaFBS/2);
    yPME1b(filt) = YpBar1(filt)+pont0*sin(pi/2/p+deltaFBS/2);
    [a2,b2,c2]=retta_per_2pti(XpBar2(filt),YpBar2(filt),xxD2k(filt),yyD2k(filt));
    a1 = tan(pi/2/p+deltaFBS/2-pi/2).*filt;
    b1 = -1.*filt;
    c1 = YpBar1(filt)+pont0*sin(pi/2/p+deltaFBS/2)-a1(filt).*(XpBar1(filt)+pont0*cos(pi/2/p+deltaFBS/2));
    [xPME2b(filt),yPME2b(filt)]=intersezione_tra_rette(a1(filt),b1(filt),c1(filt),a2(filt),b2(filt),c2(filt));
end

% xPME2b(1) = XpBar2(1);
% yPME2b(1) = YpBar2(1);
% xPME1b(1) = XpBar1(1);
% yPME1b(1) = YpBar1(1);

% PM clearance
xPME1b = xPME1b+PMclearE*cos(pi/2/p-pi/2+deltaFBS/2);
yPME1b = yPME1b+PMclearE*sin(pi/2/p-pi/2+deltaFBS/2);

b1 = calc_distanza_punti([xPME1b' yPME1b'],[xxD1k' yyD1k'])';
b2 = calc_distanza_punti([xPME2b' yPME2b'],[xxD2k' yyD2k'])';
b  = min([b1;b2],[],1)-pont0/4;
% b(1) = 0;
AreaE = b.*hc;
if FBS>0
    if flagPMFBS
        AreaEMax=AreaE;
    else
        AreaEMax=geo.AreaEMax;
    end
else
    AreaEMax=AreaE;
end

if min(PMdimE)<0
    PMdimE = -PMdimE.*b;
end

PMdimE(PMdimE>b)=b(PMdimE>b);
PMdimE(PMdimE<0)=0;
PMdimE=floor(PMdimE*100)/100;

AreaE = PMdimE.*hc;

xPME2t = xPME2b+PMdimE*cos(pi/2/p+deltaFBS/2);
yPME2t = yPME2b+PMdimE*sin(pi/2/p+deltaFBS/2);
xPME1t = xPME1b+PMdimE*cos(pi/2/p+deltaFBS/2);
yPME1t = yPME1b+PMdimE*sin(pi/2/p+deltaFBS/2);


% label
xc   = zeros(1,2*nlay);
yc   = zeros(1,2*nlay);
xmag = zeros(1,2*nlay);
ymag = zeros(1,2*nlay);
xair = zeros(1,4*nlay);
yair = zeros(1,4*nlay);
Br   = zeros(1,2*nlay);

indexZone = repmat(1:1:2,[1,nlay]);

xc(indexZone==1) = (xPMC1b+xPMC1t+xPMC2b+xPMC2t)/4;
yc(indexZone==1) = (yPMC1b+yPMC1t+yPMC2b+yPMC2t)/4;
xc(indexZone==2) = (xPME1b+xPME1t+xPME2b+xPME2t)/4;
yc(indexZone==2) = (yPME1b+yPME1t+yPME2b+yPME2t)/4;

xc(indexZone==1) = xc(indexZone==1).*(PMdimC./abs(PMdimC));
xc(indexZone==2) = xc(indexZone==2).*(PMdimE./abs(PMdimE));

xmag(indexZone==1) = ones(1,nlay);
ymag(indexZone==1) = zeros(1,nlay);
xmag(indexZone==2) = cos(pi/2/p-pi/2)*ones(1,nlay);
ymag(indexZone==2) = sin(pi/2/p-pi/2)*ones(1,nlay);

Br(indexZone==1) = mat.LayerMag.Br;
Br(indexZone==2) = mat.LayerMag.Br;

indexZone = repmat(1:1:4,[1,nlay]);

xair(indexZone==1) = (XpontRadBarDx+XpontRadBarSx+XpontRadDx+XpontRadSx)/4;   % air between radial ribs and central PM
yair(indexZone==1) = (YpontRadBarDx+YpontRadBarSx+YpontRadDx+YpontRadSx)/4;
xair(indexZone==2) = (XpontSplitBarDx(2,:)+XpontSplitBarSx(2,:)+XpontSplitDx(2,:)+XpontSplitSx(2,:))/4;  % air between split rib and central PM
yair(indexZone==2) = (YpontSplitBarDx(2,:)+YpontSplitBarSx(2,:)+YpontSplitDx(2,:)+YpontSplitSx(2,:))/4;
xair(indexZone==3) = (xPME1b+xPME2b+XpBar1)/3;  %barrier corner
yair(indexZone==3) = (yPME1b+yPME2b+YpBar1)/3;
xair(indexZone==4) = (xxD1k+xxD2k+xpont)/3;  % end barrier
yair(indexZone==4) = (yyD1k+yyD2k+ypont)/3;

% filt xair
xair(indexZone==3) = xair(indexZone==3).*(PMdimE./abs(PMdimE));     % if PME is not present, corner label must be disabled
xair(indexZone==1) = xair(indexZone==1).*(PMdimC./abs(PMdimC));     % if PMC is not present, rad rib label must be disables
filtClear = ones(size(PMclearE));
filtClear(PMclearE>0) = NaN;
xair(indexZone==3) = xair(indexZone==3).*filtClear; % if PME are with clearance, corner label must be disabled
filtClear = ones(size(PMclearC));
filtClear(PMclearC>0) = NaN;
xair(indexZone==1) = xair(indexZone==1).*filtClear; % if PMC are with clearance, rad rib labels must be disabled

%% save variables
temp.AreaCMax = AreaCMax;
temp.AreaEMax = AreaEMax;
temp.AreaC    = AreaC;
temp.AreaE    = AreaE;

geo.AreaCMax  = AreaCMax;
geo.AreaEMax  = AreaEMax;
geo.AreaC     = AreaC;
geo.AreaE     = AreaE;
geo.PMdim     = [PMdimC;PMdimE];
geo.PMclear   = [PMclearC;PMclearE];

temp.xc       = xc;
temp.yc       = yc;
temp.xair     = xair;
temp.yair     = yair;
temp.xmag     = xmag;
temp.ymag     = ymag;
temp.zmag     = zeros(size(xmag));

temp.xPMC1b = xPMC1b;
temp.yPMC1b = yPMC1b;
temp.xPMC1t = xPMC1t;
temp.yPMC1t = yPMC1t;
temp.xPMC2b = xPMC2b;
temp.yPMC2b = yPMC2b;
temp.xPMC2t = xPMC2t;
temp.yPMC2t = yPMC2t;
temp.xPME1b = xPME1b;
temp.yPME1b = yPME1b;
temp.xPME1t = xPME1t;
temp.yPME1t = yPME1t;
temp.xPME2b = xPME2b;
temp.yPME2b = yPME2b;
temp.xPME2t = xPME2t;
temp.yPME2t = yPME2t;

mat.LayerMag.Br = [Br Br];








