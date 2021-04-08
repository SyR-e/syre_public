% Copyright 2018
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

function [geo,mat,temp]=nodes_rotor_Seg(geo,mat)
%% Input
r      = geo.r;                 % Raggio del rotore al traferro
x0     = geo.x0;                % Centro fittizio
Ar     = geo.Ar;
pont0  = geo.pont0;             % minimum mechanical tolerance
pontT  = geo.pontT;             % Airgap ribs [mm]
p      = geo.p;                 % Paia poli
nlay   = geo.nlay;              % N° layers
dalpha = geo.dalpha;            % Angoli dalpha
alpha  = cumsum(dalpha); 
dx     = geo.dx;
dxIB   = geo.dxIB;
kOB    = geo.kOB;
delta_FBS   = geo.delta_FBS;    % flux barrier shift angle
hfe_min     = geo.hfe_min;
RotorFillet = geo.RotorFillet;
hcShrink = geo.hcShrink;

%% Initialitation

XcRibTraf1 = zeros(1,nlay);
XcRibTraf2 = zeros(1,nlay);
YcRibTraf1 = zeros(1,nlay);
YcRibTraf2 = zeros(1,nlay);

XpBar1 = zeros(1,nlay);
YpBar1 = zeros(1,nlay);
XpBar2 = zeros(1,nlay);
YpBar2 = zeros(1,nlay);

xxD1k = zeros(1,nlay); 
yyD1k = zeros(1,nlay);
xxD2k = zeros(1,nlay);
yyD2k = zeros(1,nlay);

XcRacc_B1  = NaN(1,nlay);
YcRacc_B1  = NaN(1,nlay);
xRaccR1_B1 = NaN(1,nlay);
yRaccR1_B1 = NaN(1,nlay);
xRaccR2_B1 = NaN(1,nlay);
yRaccR2_B1 = NaN(1,nlay);
XcRacc_B2  = NaN(1,nlay);
YcRacc_B2  = NaN(1,nlay);
xRaccR1_B2 = NaN(1,nlay);
yRaccR1_B2 = NaN(1,nlay);
xRaccR2_B2 = NaN(1,nlay);
yRaccR2_B2 = NaN(1,nlay);

%%%Tangential Rotor Fillet 
temp.xC1k  = nan(1,nlay);
temp.yC1k  = nan(1,nlay);
temp.xC2k  = nan(1,nlay);
temp.yC2k  = nan(1,nlay);
temp.xC3k  = nan(1,nlay);
temp.yC3k  = nan(1,nlay);
temp.xC4k  = nan(1,nlay);
temp.yC4k  = nan(1,nlay);
temp.xC01k = nan(1,nlay);
temp.yC01k = nan(1,nlay);
temp.xC02k = nan(1,nlay);
temp.yC02k = nan(1,nlay);

% Check array
LowDimBarrier = zeros(1,nlay);

%% Computation of the circumferences centered in (x0,0) to draw the barriers 

% Find xpont ypont
beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,r,x0);      
rbeta = (x0 - r * cos(alpha*pi/180))./(cos(beta*pi/180)); 
[xpont,ypont] = calc_intersezione_cerchi(r-pontT, rbeta, x0);

% Check
ii = (imag(xpont)~=0 | imag(ypont)~=0);
xpont(ii) = (r-2*pontT(ii)).*cos(alpha(ii)*pi/180);
ypont(ii) = (r-2*pontT(ii)).*sin(alpha(ii)*pi/180);
LowDimBarrier(ii) = 1;

Bx0 = x0-(rbeta);
geo.Bx0 = Bx0;

%% Determination of air thickness and check the feasibility of the geometry

if delta_FBS==0
    geo = calcHcCheckGeoControlwDx(geo);
    B1k = geo.B1k;
    B2k = geo.B2k;
else
    hc = geo.hc;
    B1k = Bx0-hc/2+dx.*hc/2;
    B2k = Bx0+hc/2+dx.*hc/2;
end

hc = geo.hc;

ii = find(abs(hc)<pont0/2);
B2k(ii) = min([B1k(ii),B2k(ii)])+pont0/2;
hc(ii) = B2k(ii)-B1k(ii);

%% Find Traf1 Traf2 points
m = tan(pi/2/p+delta_FBS/2)*ones(1,nlay);

x2 = xpont+(B2k-Bx0).*sin(pi/2/p+delta_FBS/2);
y2 = ypont-(B2k-Bx0).*cos(pi/2/p+delta_FBS/2);

q2 = y2-m.*x2;
[xTraf2,yTraf2] = intersezione_retta_circonferenza(0,0,(r-pontT),m,q2);

x1 = xpont+(B1k-Bx0).*sin(pi/2/p+delta_FBS/2);
y1 = ypont-(B1k-Bx0).*cos(pi/2/p+delta_FBS/2);
q1 = y1-m.*x1;
[xTraf1,yTraf1] = intersezione_retta_circonferenza(0,0,(r-pontT),m,q1);


ii = (imag(xTraf1)~=0 | imag(yTraf1)~=0 | imag(xTraf2)~=0 | imag(yTraf2)~=0 |hc<0);
xpont(ii) = (r-2*pontT(ii)).*cos(alpha(ii)*pi/180);
ypont(ii) = (r-2*pontT(ii)).*sin(alpha(ii)*pi/180);
LowDimBarrier(ii) = 1;

%% Outer branch narrowing factor [p.u.]

% Limit kOB
kOB(kOB>1) = ones(kOB>1);
kOB(kOB<0.3) = ones(kOB<0.3)*0.3;
kOB_val = (ones(1,nlay)-kOB).*hc/2;

[a,b,c] = retta_per_2pti(xTraf1,yTraf1,xTraf2,yTraf2);
[m,q,~] = retta_abc2mq(a,b,c);

a = (1+m.*m);
b = -2*xTraf1+2*m.*(q-yTraf1);
c = -kOB_val.*kOB_val+xTraf1.^2+(q-yTraf1).^2;
xTraf1 = (-b+sqrt(round((b.^2-4*a.*c),3)))./(2*a);
yTraf1 = m.*xTraf1+q;

a = (1+m.*m);
b = -2*xTraf2+2*m.*(q-yTraf2);
c = -kOB_val.*kOB_val+xTraf2.^2+(q-yTraf2).^2;
xTraf2 = (-b-sqrt(round((b.^2-4*a.*c),3)))./(2*a);
yTraf2 = m.*xTraf2+q;

%% Inner branch radial shift [mm]
hfe_min_check = hfe_min*3;

if any(dxIB)
    B1k_old=B1k;
    
    B1k=B1k+dxIB;
    B2k=B2k+dxIB;
    
    % Check first barrier
    if B2k(1)>xTraf2(1)
        B2k(1)=xTraf2(1)-pont0/2;
        B1k(1)=B2k(1)-hc(1);
    end
    
    % Check interference between last barrier and the shaft
    if B1k(end)<Ar
        B1k(end) = Ar+hfe_min_check;
        B2k(end) = B1k(end)+hc(end);
    end
    
    % Check if the barrier can fit between B1k and xTraf2
    if B1k(end)+sum(hc)+hfe_min*(nlay-1)>xTraf2(1)
        B1k(end)= xTraf2(1) - (sum(hc)+hfe_min*(nlay));
        B2k(end)= B1k(end)+hc(end);
    end
    
    % Apply the chosen dxIB
    for ii=1:(nlay-1)
        if B1k(ii)<B2k(ii+1)+hfe_min_check
            B1k(ii)=B2k(ii+1)+hfe_min_check;
            B2k(ii)=B1k(ii)+hc(ii);
            if B1k(ii)+sum(hc(1:ii))+hfe_min_check*(ii-1)>xTraf2(1)
                B1k(ii)= xTraf2(1) - (sum(hc(1:ii))+hfe_min_check*(ii-1))-pont0/2;
                B2k(ii)= B1k(ii)+hc(ii);
                B2k(ii+1)=B1k(ii)-hfe_min_check;
                B1k(ii+1)=B2k(ii+1)-hc(ii+1);
            end
            for jj=ii:-1:2          
                if B2k(jj)>B1k(jj-1)-hfe_min_check
                    B1k(jj-1)=B2k(jj)+hfe_min_check;
                    B2k(jj-1)=B1k(jj-1)+hc(jj-1);
                end
            end
        end
    end
end    

% Check if the YBar has a negative value
m = tan(pi/2/p+delta_FBS/2)*ones(1,nlay);
q = yTraf2 - m.*xTraf2;
XpBar2_lim = -q./m;
ii = B2k<XpBar2_lim;
B2k(ii) = XpBar2_lim(ii);
B1k(ii) = B2k(ii)-hc(ii);

% Check if xpont<B2k
ii = (xpont<B2k);
B2k(ii) = xpont(ii);
B1k(ii) = B2k(ii)-hc(ii);

if any(dxIB)
    dxIB = round(B1k-B1k_old,2);
end

%% Elbow definition
ii = not(LowDimBarrier);

% Compute pBar2
XpBar2(LowDimBarrier>0) = B2k(LowDimBarrier>0);
YpBar2(LowDimBarrier>0) = 0;
a1 = ones(1,nlay);
b1 = zeros(1,nlay);
c1 = -B2k;
m2 = tan(pi/2/p+delta_FBS/2)*ones(1,nlay);
a2 = m2;
b2 = -ones(1,nlay);
c2 = (yTraf2-m2.*xTraf2);
[XpBar2(ii),YpBar2(ii)] = intersezione_tra_rette(a1(ii),b1(ii),c1(ii),a2(ii),b2(ii),c2(ii));

% Shrink
hcShrink_val = hcShrink.*YpBar2;
hcShrink_val(hcShrink_val<0) = 0;
YpBar2_old = YpBar2;
YpBar2(ii) = YpBar2(ii) - hcShrink_val(ii);

flag_V = zeros(1,nlay);
flag_V(YpBar2<=geo.pontR+pont0 & hcShrink_val~=0) = 1;
YpBar2(flag_V==1) = 0;
[a,b,c] = retta_per_2pti(XpBar2,YpBar2,xTraf2,yTraf2);
[m,q,~] = retta_abc2mq(a,b,c);
hcAngle = atan(m);
B1k_old = B1k;
B1k(flag_V==1) = B2k(flag_V==1) - hc(flag_V==1)./sin(hcAngle(flag_V==1));

mTraf = -m.^-1;
qTraf = yTraf2-mTraf.*xTraf2;

a = (1+mTraf.^2);
b = -2*xTraf2+2*mTraf.*(qTraf-yTraf2);
c = -hc.*hc+xTraf2.^2+(qTraf-yTraf2).^2;
xTraf1_new = (-b-sqrt(b.^2-4*a.*c))./(2*a);
yTraf1_new = mTraf.*xTraf1_new+qTraf;

xTraf1(hcShrink_val~=0) = xTraf1_new(hcShrink_val~=0);
yTraf1(hcShrink_val~=0) = yTraf1_new(hcShrink_val~=0);

% Compute pBar1
XpBar1(LowDimBarrier>0) = B1k(LowDimBarrier>0);
YpBar1(LowDimBarrier>0) = 0;
a1 = ones(1,nlay);
b1 = zeros(1,nlay);
c1 = -B1k;
a2 = m;
b2 = -ones(1,nlay);
c2 = (yTraf1-m.*xTraf1);
[XpBar1(ii),YpBar1(ii)] = intersezione_tra_rette(a1(ii),b1(ii),c1(ii),a2(ii),b2(ii),c2(ii));


%Check hc shrink interference
if any(hcShrink_val)   
for jj = nlay:-1:2
    if(m(jj)*XpBar1(jj-1)+q(jj)<YpBar1(jj-1))
        YpBar2(jj) = YpBar1(jj-1)+hfe_min;
    end
end

flag_V = zeros(1,nlay);
flag_V(YpBar2<=geo.pontR+pont0 & hcShrink_val~=0) = 1;
[a,b,c] = retta_per_2pti(XpBar2,YpBar2,xTraf2,yTraf2);
[m,q,~] = retta_abc2mq(a,b,c);
hcAngle = atan(m);
B1k(flag_V==1) = B2k(flag_V==1) - hc(flag_V==1)./sin(hcAngle(flag_V==1));
B1k(flag_V==0) = B1k_old(flag_V==0);

mTraf = -m.^-1;
qTraf = yTraf2-mTraf.*xTraf2;

a = (1+mTraf.^2);
b = -2*xTraf2+2*mTraf.*(qTraf-yTraf2);
c = -hc.*hc+xTraf2.^2+(qTraf-yTraf2).^2;
xTraf1_new = (-b-sqrt(b.^2-4*a.*c))./(2*a);
yTraf1_new = mTraf.*xTraf1_new+qTraf;

xTraf1(hcShrink_val~=0) = xTraf1_new(hcShrink_val~=0);
yTraf1(hcShrink_val~=0) = yTraf1_new(hcShrink_val~=0);

% Compute pBar1
XpBar1(LowDimBarrier>0) = B1k(LowDimBarrier>0);
YpBar1(LowDimBarrier>0) = 0;
a1 = ones(1,nlay);
b1 = zeros(1,nlay);
c1 = -B1k;
a2 = m;
b2 = -ones(1,nlay);
c2 = (yTraf1-m.*xTraf1);
[XpBar1(ii),YpBar1(ii)] = intersezione_tra_rette(a1(ii),b1(ii),c1(ii),a2(ii),b2(ii),c2(ii));
end

% Compute D1k D2k
xxD1k(LowDimBarrier>0) = B1k(LowDimBarrier>0);
yyD1k(LowDimBarrier>0) = 0;
[xxD1k(ii),yyD1k(ii),XcRibTraf1(ii),YcRibTraf1(ii),~] = tg_cir(XpBar1(ii),YpBar1(ii),xTraf1(ii),yTraf1(ii),xpont(ii),ypont(ii));  

xxD2k(LowDimBarrier>0) = B2k(LowDimBarrier>0);
yyD2k(LowDimBarrier>0) = 0;
[xxD2k(ii),yyD2k(ii),XcRibTraf2(ii),YcRibTraf2(ii),~] = tg_cir(XpBar2(ii),YpBar2(ii),xTraf2(ii),yTraf2(ii),xpont(ii),ypont(ii));  

% Check D2k
ii = (xxD2k<=XpBar2 & XpBar2>xpont & LowDimBarrier==1);
XpBar2(ii) = xpont(ii);
B2k(ii) = xpont(ii);

ii = (xxD2k<=XpBar2 & LowDimBarrier==0);
rc = sqrt((xxD2k(ii)-XcRibTraf2(ii)).^2+(yyD2k(ii)-YcRibTraf2(ii)).^2);
xxD2k(ii) = XpBar2(ii);
yyD2k(ii) = -sqrt(rc.^2-(XpBar2(ii)-XcRibTraf2(ii)).^2)+YcRibTraf2(ii);
YpBar2(ii) = yyD2k(ii);

ii = (xxD1k<=XpBar1 & LowDimBarrier==0);
rc = sqrt((xxD1k(ii)-XcRibTraf1(ii)).^2+(yyD1k(ii)-YcRibTraf1(ii)).^2);
xxD1k(ii) = XpBar1(ii);
yyD1k(ii) = sqrt(rc.^2-(XpBar1(ii)-XcRibTraf1(ii)).^2)+YcRibTraf1(ii);
YpBar1(ii) = yyD1k(ii);

if (xxD2k(1)<xpont(1) && yyD2k(1)>=ypont(1))
    YpBar2(1)=YpBar2(1)/2;
    yyD2k(1)=YpBar2(1);
end

if any(hcShrink_val)
    hcShrink = round((YpBar2_old - YpBar2)./YpBar2_old,2);
else
    hcShrink(hcShrink_val==0) = 0;
end
%% No Vagati's Finger

if isfinite(RotorFillet)
    if sum(RotorFillet>hc/2)
        disp('Rotor Fillet limited');
        RotorFillet(1,RotorFillet>hc/2) = floor(hc(RotorFillet>hc/2)/2*10^2)/10^2;
        RotorFillet(2,RotorFillet>hc/2) = floor(hc(RotorFillet>hc/2)/2*10^2)/10^2;
    end

    m1 = (YpBar1-yyD1k)./(XpBar1-xxD1k);
    if ~isfinite(m1(1,1))
        m1(1,1) = m1(1,2);
    end
    
    q1 = YpBar1-m1.*XpBar1;
    m2 = -xpont./ypont;
    m1perp = -ones(1,nlay)./m1;
    m2perp = -ones(1,nlay)./m2;
    q2 = ypont-m2.*xpont;
    q3 = yyD2k-m1.*xxD2k;
    
    %tangential rib superior arc starts in C1k and ends in C2k with center in C01k
    xC2k = (sqrt(1+m1.*m1).*RotorFillet-q1+q2)./(m1-m2);
    yC2k = m2.*xC2k+q2;
    xC1k = (sqrt(1+m2.*m2).*RotorFillet-q2+q1)./(m2-m1);
    yC1k = m1.*xC1k+q1;
    q11 = yC1k-m1perp.*xC1k;
    q22 = yC2k-m2perp.*xC2k;
    xC01k = (q22-q11)./(m1perp-m2perp);
    yC01k = m2perp.*xC01k+q22;
    
    %tangential rib inferior arc starts in C3k and ends in C4k with center in C02k
    xC3k=(sqrt(1+m1.*m1).*RotorFillet+q3-q2)./(m2-m1);
    yC3k=m2.*xC3k+q2;
    xC4k=(sqrt(1+m2.*m2).*RotorFillet-q2+q3)./(m2-m1);
    yC4k=m1.*xC4k+q3;
    q33=yC3k-m2perp.*xC3k;
    q44=yC4k-m1perp.*xC4k;
    xC02k=(q44-q33)./(m2perp-m1perp);
    yC02k=m2perp.*xC02k+q33;
    
    %constraint on the inferior arc
    index=(xC4k<XpBar2) & (yC4k<YpBar2);
    if sum(index)>0
        disp('Inferior arc limited in rotor fillet');
        xC4k(index)=XpBar2(index);
        yC4k(index)=YpBar2(index);
        RotorFilletLimited=((m2(index)-m1(index)).*xC4k(index)+q2(index)-q3(index))./(sqrt(1+m2(index).*m2(index)));
        xC3k(index)=(sqrt(1+m1(index).*m1(index)).*RotorFilletLimited+q3(index)-q2(index))./(m2(index)-m1(index));
        yC3k(index)=m2(index).*xC3k(index)+q2(index);
    end
    
    %save nodes
    temp.xC1k = xC1k;
    temp.yC1k = yC1k;
    temp.xC2k = xC2k;
    temp.yC2k = yC2k;
    temp.xC3k = xC3k;
    temp.yC3k = yC3k;
    temp.xC4k = xC4k;
    temp.yC4k = yC4k;
    temp.xC01k = xC01k;
    temp.yC01k = yC01k;
    temp.xC02k = xC02k;
    temp.yC02k = yC02k;
end

geo.hc = hc;
temp.B1k = B1k;
temp.B2k = B2k;
temp.xpont = xpont;
temp.ypont = ypont;
temp.XpBar1 = XpBar1;
temp.YpBar1 = YpBar1;
temp.XpBar2 = XpBar2;
temp.YpBar2 = YpBar2;
temp.xxD1k = xxD1k;
temp.yyD1k = yyD1k;
temp.xxD2k = xxD2k;
temp.yyD2k = yyD2k;

%Shrink
temp.hcAngle = hcAngle;
temp.flag_segV = flag_V; 

[temp,geo] = calc_ribs_rad_fun(geo,mat,temp);

if YpBar2(1)<=0
    LowDimBarrier(1)=1;
end

if (LowDimBarrier(1)==1)
    temp.xc=(xxD1k(2:end)+xxD2k(2:end))/2;
    temp.yc=(yyD1k(2:end)+yyD2k(2:end))/2;
    
else
    temp.xc=(xxD1k+xxD2k)/2;
    temp.yc=(yyD1k+yyD2k)/2;
end


%% temp 
temp.Bx0 = Bx0;

temp.XcRacc_B1 = XcRacc_B1;
temp.YcRacc_B1 = YcRacc_B1;
temp.xRaccR1_B1 = xRaccR1_B1;
temp.yRaccR1_B1 = yRaccR1_B1;
temp.xRaccR2_B1 = xRaccR2_B1;
temp.yRaccR2_B1 = yRaccR2_B1;

temp.XcRacc_B2 = XcRacc_B2;
temp.YcRacc_B2 = YcRacc_B2;
temp.xRaccR1_B2 = xRaccR1_B2;
temp.yRaccR1_B2 = yRaccR1_B2;
temp.xRaccR2_B2 = xRaccR2_B2;
temp.yRaccR2_B2 = yRaccR2_B2;

temp.xTraf1 = xTraf1;
temp.xTraf2 = xTraf2;
temp.yTraf1 = yTraf1;
temp.yTraf2 = yTraf2;

temp.XcRibTraf1 = XcRibTraf1;
temp.YcRibTraf1 = YcRibTraf1;
temp.XcRibTraf2 = XcRibTraf2;
temp.YcRibTraf2 = YcRibTraf2;


%% PM geometry
[geo,mat,temp] = PMdefinition_Seg(geo,mat,temp);

if strcmp(mat.LayerMag.MatName,'Air')
    temp.Mag=[];
end

temp.LowDimBarrier = LowDimBarrier;
hf = [r,B1k]-[B2k,Ar]; 
geo.hf = hf;

% barrier transverse dimension (for permeance evaluation)
temp1_sk = calc_distanza_punti([mean([xxD1k' xxD2k'],2) mean([yyD1k' yyD2k'],2)],[mean([XpBar1' XpBar2'],2) mean([YpBar1' YpBar2'],2)]);
temp2_sk = calc_distanza_punti([mean([B1k' B2k'],2) mean([B1k' B2k'],2)*0],[mean([XpBar1' XpBar2'],2) mean([YpBar1' YpBar2'],2)]);
sk = temp1_sk' + temp2_sk' + hc;

%% Output

geo.sk       = sk;
geo.pbk      = geo.sk ./ geo.hc;
geo.la       = sum(geo.hc)/geo.r;
geo.lfe      = sum(geo.hf)/geo.r;
geo.ly       = (geo.R - (geo.r + geo.g + geo.lt))/geo.r;
geo.B1k      = B1k;
geo.B2k      = B2k;
geo.YpBar1   = YpBar1;
geo.hc       = hc;
geo.dxIB     = dxIB;
geo.kOB      = kOB;
geo.hcShrink = hcShrink;

geo.CentBarLength = temp.YpontSplitBarSx(2,:)*2;
geo.hrheight      = temp.XpontSplitDx(2,:)-temp.XpontSplitSx(2,:);
geo.RotorFillet   = RotorFillet;