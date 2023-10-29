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
hcShrink        = geo.hcShrink;

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
pontR     = geo.pontR;
pontSplit = geo.radial_ribs_split;


PMNcE     = geo.PMNc(2,:);
PMNcC     = geo.PMNc(1,:);
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

b  = zeros(1,nlay);
b1 = zeros(1,nlay);
b2 = zeros(1,nlay);

b_start = max(YpontRadBarDx,YpontRadBarSx);
b_end   = min(YpontSplitBarDx(2,:),YpontSplitBarSx(2,:));
b = b_end - b_start;
% b=min(b1,b2);

AreaC = (hc-PMclearC).*b;
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

AreaC = PMdimC.*(hc-PMclearC);

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

deltaFBS=0;


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

% PM clearance
xPME1b = xPME1b+PMclearE.*cos(pi/2/p-pi/2+deltaFBS/2+hcAngle);
yPME1b = yPME1b+PMclearE.*sin(pi/2/p-pi/2+deltaFBS/2+hcAngle);

b1 = calc_distanza_punti([xPME1b' yPME1b'],[xxD1k' yyD1k'])';
b2 = calc_distanza_punti([xPME2b' yPME2b'],[xxD2k' yyD2k'])';
b  = min([b1;b2],[],1)-pont0/4;
b(b<3*pont0) = 0;

% b(1) = 0;
AreaE = b.*(hc-PMclearE);
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

AreaE = PMdimE.*(hc-PMclearE);

xPME2t = xPME2b+PMdimE.*cos(pi/2/p+deltaFBS/2+hcAngle);
yPME2t = yPME2b+PMdimE.*sin(pi/2/p+deltaFBS/2+hcAngle);
xPME1t = xPME1b+PMdimE.*cos(pi/2/p+deltaFBS/2+hcAngle);
yPME1t = yPME1b+PMdimE.*sin(pi/2/p+deltaFBS/2+hcAngle);


%% Segmention
%EXTERNAL
PMNcE(PMdimE==0) = 1;
PMdimE_seg = PMdimE./PMNcE;
PMNcE(PMdimE_seg<1 & PMdimE_seg>0) = 1;
PMdimE_seg(PMdimE_seg<1 & PMdimE_seg>0) = PMdimE(PMdimE_seg<1 & PMdimE_seg>0);

[m1,q1] = calc_line2points(xPME1b,yPME1b,xPME1t,yPME1t);
[m2,q2] = calc_line2points(xPME2b,yPME2b,xPME2t,yPME2t);

xPME1s = xPME1b;
yPME1s = yPME1b;
xPME2s = xPME2b;
yPME2s = yPME2b;


for jj=1:nlay
    if PMNcE(jj)>1
        for ii=2:(PMNcE(jj)+1)
            [tmpX,tmpY] = calc_pointDistanceLine(m1(jj),q1(jj),PMdimE_seg(jj),xPME1s(ii-1,jj),yPME1s(ii-1,jj));
            xPME1s(ii,jj) = tmpX(1);
            yPME1s(ii,jj) = tmpY(1);
            [tmpX,tmpY] = calc_pointDistanceLine(m2(jj),q2(jj),PMdimE_seg(jj),xPME2s(ii-1,jj),yPME2s(ii-1,jj));
            xPME2s(ii,jj) = tmpX(1);
            yPME2s(ii,jj) = tmpY(1);
        end
    else
        xPME1s(2,jj) = xPME1t(jj);
        yPME1s(2,jj) = yPME1t(jj);
        xPME2s(2,jj) = xPME2t(jj);
        yPME2s(2,jj) = yPME2t(jj);
    end

end

% figure
% figSetting
% plot(xPME2b,yPME2b,'xk','DisplayName','PME2b')
% plot(xPME1b,yPME1b,'xb','DisplayName','PME1b')
% plot(xPME2t,yPME2t,'xg','DisplayName','PME2t')
% plot(xPME1t,yPME1t,'xc','DisplayName','PME1t')
% 
% plot(xPME2s,yPME2s,'xr','DisplayName','PME1t')
% plot(xPME1s,yPME1s,'xr','DisplayName','PME1t')
% axis equal

%CENTRAL
% PMNcC = [1 1];
PMNcC(PMdimC<0.5) = 1;
PMNcC = PMNcC/2;
% indexC = pontR==0 | pontSplit==1;
% PMNcC(indexC) = floor(PMNcC(indexC)/2)*2;
% PMNcC(indexC) = 1;
PMNcC(PMNcC==0) = 1;
PMNcC(~pontSplit) = round(PMNcC(~pontSplit));

PMdimC_seg = PMdimC./PMNcC;
PMNcC(PMdimC_seg<1 & PMdimC_seg>0) = 1;
PMdimC_seg(PMdimC_seg<1 & PMdimC_seg>0) = PMdimE(PMdimC_seg<1 & PMdimC_seg>0);

xPMC1s = xPMC1b;
yPMC1s = yPMC1b;
xPMC2s = xPMC2b;
yPMC2s = yPMC2b;


for jj=1:nlay
    if PMNcC(jj)>0.5
        if rem(PMNcC(jj)*2, 2) == 0
            for ii=2:(PMNcC(jj)+1)
                xPMC1s(ii,jj) = xPMC1b(jj);
                yPMC1s(ii,jj) = yPMC1s(ii-1,jj) + PMdimC_seg(jj);
                xPMC2s(ii,jj) = xPMC2b(jj);
                yPMC2s(ii,jj) = yPMC2s(ii-1,jj) + PMdimC_seg(jj);
            end
        else
            ii=2;
            xPMC1s(ii,jj) = xPMC1b(jj);
            yPMC1s(ii,jj) = yPMC1s(ii-1,jj) + PMdimC_seg(jj)/2;
            xPMC2s(ii,jj) = xPMC2b(jj);
            yPMC2s(ii,jj) = yPMC2s(ii-1,jj) + PMdimC_seg(jj)/2;
            for ii=3:PMNcC(jj)*2
                xPMC1s(ii,jj) = xPMC1b(jj);
                yPMC1s(ii,jj) = yPMC1s(ii-1,jj) + PMdimC_seg(jj);
                xPMC2s(ii,jj) = xPMC2b(jj);
                yPMC2s(ii,jj) = yPMC2s(ii-1,jj) + PMdimC_seg(jj);
            end
        end
    else
        xPMC1s(2,jj) = xPMC1t(jj);
        yPMC1s(2,jj) = yPMC1t(jj);
        xPMC2s(2,jj) = xPMC2t(jj);
        yPMC2s(2,jj) = yPMC2t(jj);
    end
end

PMNcC1 = PMNcC;
PMNcC1(rem(PMNcC*2,2)==1) = PMNcC1(rem(PMNcC*2,2)==1)+0.5;


% figure
% figSetting
% plot(xPMC1b,yPMC1b,'xk','DisplayName','PMC1b')
% plot(xPMC2b,yPMC2b,'xb','DisplayName','PMC2b')
% plot(xPMC2t,yPMC2t,'xg','DisplayName','PMC2t')
% plot(xPMC1t,yPMC1t,'xc','DisplayName','PMC1t')
% 
% plot(xPMC2s,yPMC2s,'xr','DisplayName','PME1t')
% plot(xPMC1s,yPMC1s,'xr','DisplayName','PME1t')
% axis equal

%% Label
xc   = zeros(1,nlay+(sum(PMNcE)));
yc   = zeros(1,nlay+(sum(PMNcE)));
xmag = zeros(1,nlay+(sum(PMNcE)));
ymag = zeros(1,nlay+(sum(PMNcE)));
xair = zeros(1,4*nlay);
yair = zeros(1,4*nlay);
Br   = zeros(1,nlay+(sum(PMNcE)));

% indexZone = repmat(1:1:2,[1,nlay]);
indexZone = [];
for ii=1:nlay
    indexZone = [indexZone ones(1,PMNcC1(ii)) 2*ones(1,PMNcE(ii))];
end

% xc(indexZone==1) = (xPMC1b+xPMC1t+xPMC2b+xPMC2t)/4;
% yc(indexZone==1) = (yPMC1b+yPMC1t+yPMC2b+yPMC2t)/4;
% xc(indexZone==2) = (xPME1b+xPME1t+xPME2b+xPME2t)/4;
% yc(indexZone==2) = (yPME1b+yPME1t+yPME2b+yPME2t)/4;

% xc(:,indexZone==1) = xc(indexZone==1).*(PMdimC./abs(PMdimC)).*ones(length(xc(:,1)),1);
% xc(:,indexZone==2) = xc(indexZone==2).*(PMdimE./abs(PMdimE)).*ones(length(xc(:,1)),1);

%%%%
% tmpX = xc(indexZone==2);
% tmpY = yc(indexZone==2);
for jj=1:nlay
        for ii=1:PMNcE(jj)
            tmpXE(ii,jj) = (xPME2s(ii,jj)+xPME1s(ii,jj)+xPME2s(ii+1,jj)+xPME1s(ii+1,jj))/4;
            tmpYE(ii,jj) = (yPME2s(ii,jj)+yPME1s(ii,jj)+yPME2s(ii+1,jj)+yPME1s(ii+1,jj))/4;
            angX(ii,jj) = cos(hcAngle(jj)+pi/2/p-pi/2);
            angY(ii,jj) = sin(hcAngle(jj)+pi/2/p-pi/2);
        end

        for ii=1:PMNcC1(jj)
            tmpXC(ii,jj) = (xPMC2s(ii,jj)+xPMC1s(ii,jj)+xPMC2s(ii+1,jj)+xPMC1s(ii+1,jj))/4;
            tmpYC(ii,jj) = (yPMC2s(ii,jj)+yPMC1s(ii,jj)+yPMC2s(ii+1,jj)+yPMC1s(ii+1,jj))/4;
        end
end

% remove the central PM segment
% xPMC1s(1,PMNcC==1 | rem(PMNcC, 2) == 0) = nan;
% yPMC1s(1,PMNcC==1 | rem(PMNcC, 2) == 0) = nan;
% xPMC2s(1,PMNcC==1 | rem(PMNcC, 2) == 0) = nan;
% yPMC2s(1,PMNcC==1 | rem(PMNcC, 2) == 0) = nan;


tmpXE1 = tmpXE; 
tmpXC1 = tmpXC; 
tmpXE = reshape(tmpXE(tmpXE>0),[1,(sum((indexZone==2)))]);
tmpYE = reshape(tmpYE(tmpXE1>0),[1,(sum((indexZone==2)))]);
tmpXC = reshape(tmpXC(tmpXC1>0),[1,(sum((indexZone==1)))]);
tmpYC = reshape(tmpYC(tmpXC1>0),[1,(sum((indexZone==1)))]);
angX  = reshape(angX(tmpXE1>0),[1,(sum((indexZone==2)))]);
angY  = reshape(angY(tmpXE1>0),[1,(sum((indexZone==2)))]);
tmpXE(PMdimE==0) = NaN;
tmpXC(PMdimC==0) = NaN;

xc(indexZone==2) = tmpXE;
yc(indexZone==2) = tmpYE;
xc(indexZone==1) = tmpXC;
yc(indexZone==1) = tmpYC;
xmag(indexZone==2) = angX;
ymag(indexZone==2) = angY;
%%%%

xmag(indexZone==1) = ones(1,(sum((indexZone==1))));
ymag(indexZone==1) = zeros(1,(sum((indexZone==1))));

% xmag(indexZone==2) = cos(hcAngle+pi/2/p-pi/2);
% ymag(indexZone==2) = sin(hcAngle+pi/2/p-pi/2);

% Br(indexZone==1) = mat.LayerMag.Br;
% Br(indexZone==2) = mat.LayerMag.Br;
Br(indexZone==1) = mat.LayerMag.Br(1)*ones(1,(sum(indexZone==1)));
Br(indexZone==2) = mat.LayerMag.Br(1)*ones(1,(sum(indexZone==2)));

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

mirrorFlagAir = ones(1,length(xair));


tmp = mirrorFlagAir(indexZone==2);
tmp(PMclearC>0) = 0;
mirrorFlagAir(indexZone==2) = tmp;

tmp = mirrorFlagAir(indexZone==4);
tmp(PMclearE>0 & PMclearC>0 & pontR==0 & hcShrink==0) = 0;
% tmp(PMclearE>0 & PMclearC>0 & pontR==0) = 0;
mirrorFlagAir(indexZone==4) = tmp;

% mirrorFlagAir(isnan(xair)) = 0; 

%mirror flag PM label
zz = 1;
flagOdd = 0;
clear NcC
for ii=1:nlay
    for jj=1:PMNcC1(ii)
        NcC(zz+jj-1) = PMNcC(ii);
        if flagOdd
            NcC(zz+jj-1) = NcC(zz+jj-1)+0.5;
        end
        if rem(PMNcC(ii)*2, 2) == 1
            flagOdd = 1;
        end
        %         if (pontR(ii)>0 && pontSplit(ii)==0)
        %             NcC(zz+jj-1) = 1000;
        %         end
    end
    flagOdd = 0;
    for kk=1:PMNcE(ii)
        NcC(zz+jj+kk-1) = 0;
    end
    zz = length(NcC)+1;
end
mirrorFlag = ones(1,length(xc));
% mirrorFlag(((NcC==1 | rem(NcC, 2) == 0) & NcC>0) & NcC < 1000) = 0;
mirrorFlag((NcC*2)==1) = 0;
mirrorFlag(rem(NcC*2, 2) == 1) = 0;

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
geo.PMNc      = [PMNcC*2; PMNcE];
geo.PMclear   = [PMclearC;PMclearE];

temp.xc         = xc;
temp.yc         = yc;
temp.xair       = xair;
temp.yair       = yair;
temp.xmag       = xmag;
temp.ymag       = ymag;
temp.zmag       = zeros(size(xmag));
temp.mirrorFlag    = mirrorFlag;
temp.mirrorFlagAir = mirrorFlagAir;

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
temp.xPME2s  = xPME2s;
temp.yPME2s  = yPME2s;
temp.xPME1s  = xPME1s;
temp.yPME1s  = yPME1s;
temp.xPMC1s  = xPMC1s;
temp.yPMC1s  = yPMC1s;
temp.xPMC2s  = xPMC2s;
temp.yPMC2s  = yPMC2s;
temp.hcAngle = hcAngle;


mat.LayerMag.Br = [Br Br];