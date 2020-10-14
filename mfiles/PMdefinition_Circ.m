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

function [geo,mat,temp] = PMdefinition_Circ(geo,mat,temp)
% 
% [geo,temp] = PMdefinition_Circ(geo,temp,mat)
%

% barrier
xxD1k  = temp.xxD1k;
yyD1k  = temp.yyD1k;
xxD2k  = temp.xxD2k;
yyD2k  = temp.yyD2k;
xpont  = temp.xpont;
ypont  = temp.ypont;
B1k    = temp.B1k;
B2k    = temp.B2k;

% Rad Rib
YpontRadSx    = temp.YpontRadSx;
XpontRadSx    = temp.XpontRadSx;
YpontRadDx    = temp.YpontRadDx;
XpontRadDx    = temp.XpontRadDx;
XpontRadBarDx = temp.XpontRadBarDx;
YpontRadBarDx = temp.YpontRadBarDx;
XpontRadBarSx = temp.XpontRadBarSx;
YpontRadBarSx = temp.YpontRadBarSx;


% from geo
%p         = geo.p;
nlay      = geo.nlay;
hc        = geo.hc;
FBS       = abs(geo.delta_FBS);
%betaC     = geo.PMdim(1,:);
%betaE     = geo.PMdim(2,:);
PMdim     = geo.PMdim(1,:);
flagPMFBS = geo.flagPMFBS;
x0        = geo.x0;
PMdir     = geo.PMdir;
PMclear   = geo.PMclear(1,:);
pont0     = geo.pont0;

% initialize variables for PM area.
%   C --> central portion of the barrier
%   E --> external portion of the barrier (arms)
AreaCMax = zeros(1,nlay);
AreaEMax = zeros(1,nlay);
AreaC    = zeros(1,nlay);
AreaE    = zeros(1,nlay);


% compute the PM area (2 PM section each half barrier for PM direction)
eta1 = atan2(yyD1k,(x0-xxD1k));
eta2 = atan2(yyD2k,(x0-xxD2k));
etaMax = min([eta1;eta2],[],1);

etaRib = atan2(YpontRadBarDx,(x0-XpontRadBarDx));

AreaC = pi*(etaMax-etaRib)/(2*pi).*((x0-B1k).^2-(x0-B2k).^2);
b     = AreaC./hc;
if FBS>0
    if flagPMFBS
        AreaCMax = AreaC;
    else
        AreaCMax = geo.AreaCMax;
    end
else
    AreaCMax = AreaC;
end

% AreaC = AreaCMax.*PMdim./b;
if min(PMdim)<0 % per-unit input (during optimization)
    PMdim = -PMdim.*b;
end

PMdim(PMdim>b)=b(PMdim>b);
PMdim(PMdim<0)=0;
PMdim=floor(PMdim*100)/100;

AreaC = PMdim.*hc;


% Check clearance size: the minimum PM thickness is set to pont0
PMclearLim = hc-pont0*ones(size(hc));
PMclear(PMclear>PMclearLim) = PMclearLim(PMclear>PMclearLim);

xPMC2b = XpontRadBarDx;
yPMC2b = YpontRadBarDx;
m = -tan(etaRib);
q = 0;
[xPMC1b,yPMC1b] = intersezione_retta_circonferenza(0,0,x0-B1k-PMclear,m,q);
xPMC1b = x0-xPMC1b;
yPMC1b = -yPMC1b;

m = -tan((etaMax.*PMdim./b+etaRib)/2);
q = 0;
[xPMC1t,yPMC1t] = intersezione_retta_circonferenza(0,0,x0-B1k-PMclear,m,q);
xPMC1t = x0-xPMC1t;
yPMC1t = -yPMC1t;
[xPMC2t,yPMC2t] = intersezione_retta_circonferenza(0,0,x0-B2k,m,q);
xPMC2t = x0-xPMC2t;
yPMC2t = -yPMC2t;

xPME1b = xPMC1t;
yPME1b = yPMC1t;
xPME2b = xPMC2t;
yPME2b = yPMC2t;

m = -tan(etaMax.*PMdim./b);
q = 0;
[xPME2t,yPME2t] = intersezione_retta_circonferenza(0,0,x0-B2k,m,q);
xPME2t = x0-xPME2t;
yPME2t = -yPME2t;
[xPME1t,yPME1t] = intersezione_retta_circonferenza(0,0,x0-B1k-PMclear,m,q);
xPME1t = x0-xPME1t;
yPME1t = -yPME1t;

% label

xc   = zeros(1,2*nlay);
yc   = zeros(1,2*nlay);
xair = zeros(1,2*nlay);
yair = zeros(1,2*nlay);
Br   = zeros(1,2*nlay);

index=1:1:2*nlay;

indexZone = repmat(1:1:2,[1,nlay]);

xc(indexZone==1) = (xPMC1b+xPMC1t+xPMC2b+xPMC2t)/4;
yc(indexZone==1) = (yPMC1b+yPMC1t+yPMC2b+yPMC2t)/4;
xc(indexZone==2) = (xPME1b+xPME1t+xPME2b+xPME2t)/4;
yc(indexZone==2) = (yPME1b+yPME1t+yPME2b+yPME2t)/4;

xc(indexZone==1) = xc(indexZone==1).*(PMdim./abs(PMdim));
xc(indexZone==2) = xc(indexZone==2).*(PMdim./abs(PMdim));

xair(indexZone==1) = (XpontRadBarDx+XpontRadBarSx+XpontRadDx+XpontRadSx)/4;
yair(indexZone==1) = (YpontRadBarDx+YpontRadBarSx+YpontRadDx+YpontRadSx)/4;
xair(indexZone==2) = (xxD1k+xxD2k+xpont)/3;
yair(indexZone==2) = (yyD1k+yyD2k+ypont)/3;

xair(indexZone==1) = xair(indexZone==1).*(PMdim./abs(PMdim));
filtClear = ones(size(PMclear));
filtClear(PMclear>0) = NaN;
xair(indexZone==1) = xair(indexZone==1).*filtClear; % if PMC are with clearance, rad rib labels must be disabled

if PMdir=='p'   % PM magnetization parallel to q-axis
    xmag = ones(1,2*nlay);
    ymag = zeros(1,2*nlay);
else
    etaMag = atan2(yc,xc-x0);
    xmag = -cos(etaMag);
    ymag = -sin(etaMag);
end

zmag = zeros(1,2*nlay);

Br(rem(index,2)==1) = mat.LayerMag.Br;
Br(rem(index,2)==0) = mat.LayerMag.Br;

%% save variables
temp.AreaCMax = AreaCMax;
temp.AreaEMax = AreaEMax;
temp.AreaC    = AreaC;
temp.AreaE    = AreaE;

geo.AreaCMax  = AreaCMax;
geo.AreaEMax  = AreaEMax;
geo.AreaC     = AreaC;
geo.AreaE     = AreaE;
geo.PMdim     = [PMdim;zeros(size(PMdim))];
geo.PMclear   = [PMclear;zeros(size(PMclear))];

temp.xc       = xc;
temp.yc       = yc;
temp.xair     = xair;
temp.yair     = yair;
temp.xmag     = xmag;
temp.ymag     = ymag;
temp.zmag     = zmag;

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

