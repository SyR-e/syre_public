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
XpontRadBarDx = abs(temp.XpontRadBarDx);
YpontRadBarDx = abs(temp.YpontRadBarDx);
XpontRadBarSx = abs(temp.XpontRadBarSx);
YpontRadBarSx = abs(temp.YpontRadBarSx);

% Split Rad Rib
XpontSplitBarSx = abs(temp.XpontSplitBarSx);
YpontSplitBarSx = abs(temp.YpontSplitBarSx);
XpontSplitBarDx = abs(temp.XpontSplitBarDx);
YpontSplitBarDx = abs(temp.YpontSplitBarDx);
XpontSplitDx    = temp.XpontSplitDx;
YpontSplitDx    = temp.YpontSplitDx;
XpontSplitSx    = temp.XpontSplitSx;
YpontSplitSx    = temp.YpontSplitSx;

hcAngle         = temp.hcAngle-pi/2/geo.p;

% flag_cent= temp.flag_cent;
% flag_ext= temp.flag_ext;
% flag_shiftUP= temp.flag_shiftUP;

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

% if ~isempty(geo.RQnames)
%     PMdimC(2*PMdimC<hc & geo.radial_ribs_split) = 0;
%     PMdimC(PMdimC<hc & ~geo.radial_ribs_split) = 0;
%     PMdimE(PMdimE<hc) = 0;
% end

%ang_rib   = geo.ang_rib;

% initialize variables for PM area.
%   C --> central portion of the barrier
%   E --> external portion of the barrier (arms)

% initialize points for PM definition
%   1 --> inner side of the barrier (1==sx)
%   2 --> outer side of the barrier (2==dx)
%   t --> top side of the magnet
%   b --> bottom side of the magnet

% point definition - central section

% b = YpontSplitBarDx(2,:)-YpontRadBarDx-pont0/4;
% flagDx = YpontSplitBarDx(2,:)<YpontSplitBarSx(2,:);

b = zeros(1,nlay);
b1 = zeros(1,nlay);
b2 = zeros(1,nlay);
% b(flagDx==1) = YpontSplitBarDx(2,flagDx==1)-YpontRadBarDx(flagDx==1)-pont0/4;
% b(flagDx==0) = YpontSplitBarSx(2,flagDx==0)-YpontRadBarSx(flagDx==0)-pont0/4;

% b1 = YpontSplitBarDx(2,:)-YpontRadBarDx-pont0/4;
% b2 = YpontSplitBarSx(2,:)-YpontRadBarSx-pont0/4;

b_start = max(YpontRadBarDx,YpontRadBarSx);
b_end   = min(YpontSplitBarDx(2,:),YpontSplitBarSx(2,:));
b = b_end - b_start;
% b=min(b1,b2);

%if pont!!
% index_PM1= flag_cent==1 & flag_shiftUP==0;
% if sum(index_PM1)>0
%     index_ang1= ang_rib(index_PM1)<0;
%     if sum(index_ang1)>0
%         b(index_ang1)=  YpontSplitBarDx(2,(index_ang1))-YpontRadBarDx(index_ang1)-pont0/4;
%     end
%     index_ang1=not(index_ang1);
%     if sum(index_ang1)>0
%         b(index_ang1)=  YpontSplitBarSx(2,(index_ang1))-YpontRadBarSx(index_ang1)-pont0/4;
%     end
% end

% index_PM2= flag_shiftUP==1;
% if sum(index_PM2)>0
%         b(index_PM2)=  YpBar2(index_PM2)-YpontRadBarDx(index_PM2)-pont0/4;
% end


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
% PMdimC=floor(PMdimC*100)/100;

AreaC = PMdimC.*hc;


% Check clearance size: the minimum PM thickness is set to pont0
PMclearLim = hc-pont0*ones(size(hc));
PMclearC(PMclearC>PMclearLim) = PMclearLim(PMclearC>PMclearLim);
PMclearE(PMclearE>PMclearLim) = PMclearLim(PMclearE>PMclearLim);

xPMC1b = XpontRadBarSx+PMclearC;
yPMC1b(YpontRadBarSx<=YpontRadBarDx) = YpontRadBarDx(YpontRadBarSx<=YpontRadBarDx);
yPMC1b(YpontRadBarSx>YpontRadBarDx) = YpontRadBarSx(YpontRadBarSx>YpontRadBarDx);
xPMC2b = XpontRadBarDx;
yPMC2b(YpontRadBarSx<=YpontRadBarDx) = YpontRadBarDx(YpontRadBarSx<=YpontRadBarDx);
yPMC2b(YpontRadBarSx>YpontRadBarDx) = YpontRadBarSx(YpontRadBarSx>YpontRadBarDx);

xPMC2t = xPMC2b;
yPMC2t = yPMC2b+PMdimC;
xPMC1t = xPMC1b;
yPMC1t = yPMC2t;

% point definition - external section
xPME2b = XpontSplitBarDx(1,:)+pont0*cos(pi/2/p+deltaFBS/2+hcAngle);
yPME2b = YpontSplitBarDx(1,:)+pont0*sin(pi/2/p+deltaFBS/2+hcAngle);

xPME2b(yPME2b<YpBar2)=XpBar2(yPME2b<YpBar2)+pont0*cos(pi/2/p+deltaFBS/2+hcAngle(yPME2b<YpBar2));
yPME2b(yPME2b<YpBar2)=YpBar2(yPME2b<YpBar2)+pont0*sin(pi/2/p+deltaFBS/2+hcAngle(yPME2b<YpBar2));


[a1,b1,c1] = retta_per_2pti(XpBar1,YpBar1,xxD1k,yyD1k);
a2 = tan(pi/2/p+deltaFBS/2+hcAngle-pi/2).*ones(1,nlay);
b2 = -1*ones(1,nlay);
%c2 = YpontSplitBarDx(1,:)+pont0*sin(pi/2/p+deltaFBS/2)-a2.*(XpontSplitBarDx(1,:)+pont0*cos(pi/2/p+deltaFBS/2));
c2 = yPME2b -a2.*xPME2b;
[xPME1b,yPME1b]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);

filt=YpBar1>yPME1b;
if sum(filt)
    xPME1b(filt) = XpBar1(filt)+pont0*cos(pi/2/p+deltaFBS/2+hcAngle(filt));
    yPME1b(filt) = YpBar1(filt)+pont0*sin(pi/2/p+deltaFBS/2+hcAngle(filt));
%   [a2,b2,c2]=retta_per_2pti(XpBar2(filt),YpBar2(filt),xxD2k(filt),yyD2k(filt));
    [a2,b2,c2]=retta_per_2pti(XpBar2,YpBar2,xxD2k,yyD2k);
    a1 = tan(pi/2/p+deltaFBS/2+hcAngle-pi/2).*filt;
    b1 = -1.*filt;
    c1(filt) = (YpBar1(filt)+pont0*sin(pi/2/p+deltaFBS/2+hcAngle(filt))-a1(filt).*(XpBar1(filt)+pont0*cos(pi/2/p+deltaFBS/2+hcAngle(filt))));
    [xPME2b(filt),yPME2b(filt)]=intersezione_tra_rette(a1(filt),b1(filt),c1(filt),a2(filt),b2(filt),c2(filt));
end

filt=YpontSplitBarSx(1,:)>yPME1b;
if sum(filt)
    xPME1b(filt) = XpontSplitBarSx(1,filt)+pont0*cos(pi/2/p+deltaFBS/2+hcAngle(filt));
    yPME1b(filt) = YpontSplitBarSx(1,filt)+pont0*sin(pi/2/p+deltaFBS/2+hcAngle(filt));
%     [a2,b2,c2]=retta_per_2pti(XpBar2(filt),YpBar2(filt),xxD2k(filt),yyD2k(filt));
    [a2,b2,c2]=retta_per_2pti(XpontSplitBarDx(1,filt),YpontSplitBarDx(1,filt),xxD2k,yyD2k);
    [a3,b3,c3]=retta_per_2pti(XpBar2,YpBar2,xxD2k,yyD2k);
    a1 = tan(pi/2/p+deltaFBS/2+hcAngle-pi/2).*filt;
    b1 = -1.*filt;
    c1(filt) = YpontSplitBarSx(1,filt)+pont0*sin(pi/2/p+deltaFBS/2+hcAngle(filt))-a1(filt).*(XpontSplitBarSx(1,filt)+pont0*cos(pi/2/p+deltaFBS/2+hcAngle(filt)));
    [xPME2b(filt),yPME2b(filt)]=intersezione_tra_rette(a1(filt),b1(filt),c1(filt),a2(filt),b2(filt),c2(filt));
end



% xPME2b(1) = XpBar2(1);
% yPME2b(1) = YpBar2(1);
% xPME1b(1) = XpBar1(1);
% yPME1b(1) = YpBar1(1);


% xPME2t(yPME2b<YpBar2)=XpBar1(yPME2b<YpBar2);
% yPME2t(yPME2b<YpBar2)=YpBar2(yPME2b<YpBar2);

% 
% xPME2t = xPME2b+PMdimE*cos(pi/2/p+deltaFBS/2);
% yPME2t = yPME2b+PMdimE*sin(pi/2/p+deltaFBS/2);


% XpontSplitBarDxRot(1,:)=XpontSplitBarDx(2,:)*cos(pi/2/p+deltaFBS/2)-YpontSplitBarDx(2,:)*sin(pi/2/p+deltaFBS/2);
% YpontSplitBarDxRot(1,:)=XpontSplitBarDx(2,:)*sin(pi/2/p+deltaFBS/2)+YpontSplitBarDx(2,:)*cos(pi/2/p+deltaFBS/2);
% 
% XpontSplitBarSxRot(1,:)=XpontSplitBarSx(2,:)*cos(pi/2/p+deltaFBS/2)-YpontSplitBarSx(2,:)*sin(pi/2/p+deltaFBS/2);
% YpontSplitBarSxRot(1,:)=XpontSplitBarSx(2,:)*sin(pi/2/p+deltaFBS/2)+YpontSplitBarSx(2,:)*cos(pi/2/p+deltaFBS/2);
% 
% flagDx = zeros(1,nlay);
% flagDx = YpontSplitBarDxRot(1,:)<YpontSplitBarSxRot(1,:);
% 
% limPMext = zeros(1,nlay);
% XlimPMext(flagDx==1) = XpontSplitBarDx(2,flagDx==1);
% YlimPMext(flagDx==1) = YpontSplitBarDx(2,flagDx==1);
% XlimPMext(flagDx==0) = XpontSplitBarSx(2,flagDx==0);
% YlimPMext(flagDx==0) = YpontSplitBarSx(2,flagDx==0);
% 
% index_PMlimext= YlimPMext< YpBar2;
% if sum(index_PMlimext)
%        XlimPMext(index_PMlimext)=XpBar2(index_PMlimext);
%        YlimPMext(index_PMlimext)=YpBar2(index_PMlimext);
% end 
% 
% b=zeros(1,nlay);
% b(flagDx==1) = calc_distanza_punti([XlimPMext(flagDx==1)' YlimPMext(flagDx==1)'],[xxD2k(flagDx==1)' yyD2k(flagDx==1)'])';
% b(flagDx==0) = calc_distanza_punti([XlimPMext(flagDx==0)' YlimPMext(flagDx==0)'],[xxD1k(flagDx==0)' yyD1k(flagDx==0)'])';
    
% index_PMe1= flag_ext==1 & ang_rib<0
% if sum(index_PMe1)
%     xPME2b(index_PMe1) = XpontSplitBarDx(1,index_PMe1)+pont0*cos(pi/2/p+deltaFBS/2);
%     yPME2b(index_PMe1) = YpontSplitBarDx(1,index_PMe1)+pont0*sin(pi/2/p+deltaFBS/2);
%     [a1,b1,c1] = retta_per_2pti(XpBar1,YpBar1,xxD1k,yyD1k);
%     a2 = tan(pi/2/p+deltaFBS/2-pi/2).*(index_PMe1);
%     b2 = -1*(index_PMe1);
%     c2 = YpontSplitBarDx(1,(index_PMe1))+pont0*sin(pi/2/p+deltaFBS/2)-a2(index_PMe1).*(XpontSplitBarDx(1,(index_PMe1))+pont0*cos(pi/2/p+deltaFBS/2)).*(index_PMe1);
%     [xPME1b(index_PMe1),yPME1b(index_PMe1)]=intersezione_tra_rette(a1(index_PMe1),b1(index_PMe1),c1(index_PMe1),a2(index_PMe1),b2(index_PMe1),c2(index_PMe1));
% end
% 
% index_PMe1= flag_ext==1 & ang_rib>0;
% if sum(index_PMe1)
%         xPME1b(index_PMe1) = XpontSplitBarSx(1,index_PMe1)+pont0*cos(pi/2/p+deltaFBS/2);
%         yPME1b(index_PMe1) = YpontSplitBarSx(1,index_PMe1)+pont0*sin(pi/2/p+deltaFBS/2);
%         
%     [a1,b1,c1] = retta_per_2pti(XpBar2,YpBar2,xxD2k,yyD2k);
%     a2 = tan(pi/2/p+deltaFBS/2-pi/2).*(index_PMe1);
%     b2 = -1*(index_PMe1);
%     c2 = YpontSplitBarSx(1,(index_PMe1))+pont0*sin(pi/2/p+deltaFBS/2)-a2(index_PMe1).*(XpontSplitBarSx(1,(index_PMe1))+pont0*cos(pi/2/p+deltaFBS/2)).*(index_PMe1);
%     [xPME2b(index_PMe1),yPME2b(index_PMe1)]=intersezione_tra_rette(a1(index_PMe1),b1(index_PMe1),c1(index_PMe1),a2(index_PMe1),b2(index_PMe1),c2(index_PMe1));
% end



% xPME2b(1) = XpBar2(1);
% yPME2b(1) = YpBar2(1);
% xPME1b(1) = XpBar1(1);
% yPME1b(1) = YpBar1(1);

% PM clearance
xPME1b = xPME1b+PMclearE.*cos(pi/2/p-pi/2+deltaFBS/2+hcAngle);
yPME1b = yPME1b+PMclearE.*sin(pi/2/p-pi/2+deltaFBS/2+hcAngle);

b1 = calc_distanza_punti([xPME1b' yPME1b'],[xxD1k' yyD1k'])';
b2 = calc_distanza_punti([xPME2b' yPME2b'],[xxD2k' yyD2k'])';
b  = min([b1;b2],[],1)-pont0/4;
b(b<3*pont0) = 0;

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

if ~isempty(geo.RQnames)
    PMdimC(2*PMdimC<hc & geo.radial_ribs_split) = 0;
    PMdimC(PMdimC<hc & ~geo.radial_ribs_split) = 0;
    PMdimE(PMdimE<hc) = 0;
end

PMdimE(PMdimE>b)=b(PMdimE>b);
PMdimE(PMdimE<0)=0;
PMdimE=floor(PMdimE*100)/100;

AreaE = PMdimE.*hc;

xPME2t = xPME2b+PMdimE.*cos(pi/2/p+deltaFBS/2+hcAngle);
yPME2t = yPME2b+PMdimE.*sin(pi/2/p+deltaFBS/2+hcAngle);
xPME1t = xPME1b+PMdimE.*cos(pi/2/p+deltaFBS/2+hcAngle);
yPME1t = yPME1b+PMdimE.*sin(pi/2/p+deltaFBS/2+hcAngle);


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
% xmag(indexZone==2) = cos(pi/2/p-pi/2)*ones(1,nlay);
% ymag(indexZone==2) = sin(pi/2/p-pi/2)*ones(1,nlay);
xmag(indexZone==2) = cos(hcAngle+pi/2/p-pi/2);
ymag(indexZone==2) = sin(hcAngle+pi/2/p-pi/2);

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

temp.xPMC1b  = xPMC1b;
temp.yPMC1b  = yPMC1b;
temp.xPMC1t  = xPMC1t;
temp.yPMC1t  = yPMC1t;
temp.xPMC2b  = xPMC2b;
temp.yPMC2b  = yPMC2b;
temp.xPMC2t  = xPMC2t;
temp.yPMC2t  = yPMC2t;
temp.xPME1b  = xPME1b;
temp.yPME1b  = yPME1b;
temp.xPME1t  = xPME1t;
temp.yPME1t  = yPME1t;
temp.xPME2b  = xPME2b;
temp.yPME2b  = yPME2b;
temp.xPME2t  = xPME2t;
temp.yPME2t  = yPME2t;
temp.hcAngle = hcAngle;

mat.LayerMag.Br = [Br Br];
