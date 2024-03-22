% Copyright 2019
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

function [geo,mat,temp]=nodes_rotor_Vtype_v3(geo,mat)

%% Input data
r         = geo.r;                  % rotor outer radius [mm]
Ar        = geo.Ar;                 % shaft radius [mm]
l         = geo.l;                  % stack length [mm]
g         = geo.g;                  % airgap length [mm]
pont0     = geo.pont0;              % mechanical tolerance [mm]
pontT     = geo.pontT;              % tangential ribs thickness [mm]
pontR     = geo.pontR;              % radial ribs thinkness [mm]
beta_pu   = geo.betaPMshape;        % per-unit PM angle
p         = geo.p;                  % pole pairs number
nlay      = geo.nlay;               % number of rotor layers
dalpha_pu = geo.dalpha_pu;          % rotor slot positions [pu]
hc_pu     = geo.hc_pu;              % barrier thickness [pu]
hfe_min   = geo.hfe_min;            % minimum iron thickness between two barriers [mm]
dx        = geo.dx;                 % barrier end shape factor
nmax      = geo.nmax;               % maximum speed [rpm]
PMdimC    = geo.PMdim(1,:);         % PM dimensions [mm]
PMclearC  = geo.PMclear(1,:);       % PM clearance [mm]
rhoFE     = mat.Rotor.kgm3;         % Densità del ferro di rotore [kg/m3]
sigmaFe   = mat.Rotor.sigma_max;    % Yield strength of the rotor iron [MPa]

radial_ribs_eval = geo.radial_ribs_eval;    % radial ribs flag
RotorFillet1 = geo.RotorFillet1;
RotorFillet2 = geo.RotorFillet2;
% RotorFillet=geo.RotorFillet;

%% Initializations
% top barrier points (near airgap)
XcRibTraf1 = zeros(1,nlay);
XcRibTraf2 = zeros(1,nlay);
YcRibTraf1 = zeros(1,nlay);
YcRibTraf2 = zeros(1,nlay);
xxD1k      = zeros(1,nlay);
yyD1k      = zeros(1,nlay);
xxD2k      = zeros(1,nlay);
yyD2k      = zeros(1,nlay);
% bottom barrier (d-axis)
xxB1k = zeros(1,nlay);
yyB1k = zeros(1,nlay);
xxB2k = zeros(1,nlay);
yyB2k = zeros(1,nlay);

%% Check hc_pu
% min and max

hcPuMin = pont0/g;
hcPuMax = (r-Ar)/g;

hc_pu(hc_pu>hcPuMax) = hcPuMax;

hcPuLim = 20*hcPuMax;
if hcPuLim>((r-Ar)/g)
    hcPuLim=((r-Ar)/g);
end
if sum(hc_pu)>hcPuLim
    hc_pu = hc_pu/sum(hc_pu)*hcPuLim;
end

hc_pu(hc_pu<hcPuMin) = hcPuMin;

alpha = cumsum(dalpha_pu*pi/2/p);
beta  = pi/2-beta_pu*(pi/2-pi/2/p);
hc    = hc_pu*g;

xpont = (r-pontT).*cos(alpha);
ypont = (r-pontT).*sin(alpha);

% hc for more than one barrier
alpha_pont=atan2(ypont,xpont);
xcbar = (r-pontT-hc/2).*cos(alpha_pont);
ycbar = (r-pontT-hc/2).*sin(alpha_pont);

for ii=1:nlay-1
   d = ((xcbar(ii+1)-xcbar(ii)).^2+(ycbar(ii+1)-ycbar(ii)).^2).^0.5;
   if d<(hc(ii)/2+hc(ii+1)/2+hfe_min)
       dhc = (-hc(ii)/2-hc(ii+1)/2-3*hfe_min+d);			%MODIFICA
%        dhc = (hc(ii)/2-hc(ii+1)/2+hfe_min-d)/2;
       hc(ii:ii+1) = hc(ii:ii+1)+dhc.*hc(ii:ii+1)/sum(hc(ii:ii+1));
       disp('#1 space between barrier end')
   end
end

hc_pu = hc/g;
geo.hc = hc;
geo.hc_pu = hc_pu;

% hc <vs> spider
alpha_pont=atan2(ypont,xpont);
xcbar = (r-pontT-hc/2).*cos(alpha_pont);
ycbar = (r-pontT-hc/2).*sin(alpha_pont);
% retta 1 = limite polo, retta 2 = retta perpendicolare passante per (xc,yc);  
m1 = tan(pi/2/p);
q1 = 0;
m2 = -1/m1;
q2 = ycbar-m2*xcbar;
[xSpider,ySpider] = intersezione_tra_rette(m1,-1,q1,m2,-1,q2);
% dSpider = ((xcbar-xSpider)^2+(ycbar-ySpider)^2)^0.5;
% xhc = xcbar-abs((hc/2)*cos(atan(m2)));
% yhc = xhc*m2+q2;
xS = xSpider+abs(hfe_min*cos(atan(m2)));
yS = xS*m2+q2;
hcMaxSpider = 2*((xcbar-xS).^2+(ycbar-yS).^2).^0.5;

if hc(nlay)>hcMaxSpider(nlay)
    if hcMaxSpider(nlay)>(hcPuMin*g)
        disp('#1 thin spider, hc limit')
        hc(nlay)=hcMaxSpider(nlay);
        %hc_pu=hc/(2*hc_half_max);
        hc_pu = hc/g;
        xcbar = (r-pontT-hc/2).*cos(alpha_pont);
        ycbar = (r-pontT-hc/2).*sin(alpha_pont);
    else
        disp('#1 thin spider, alpha limit')
        hc(nlay)=hcPuMin*g;
        hc_pu = hc/g;
        alpha_pont(end) = pi/2/p-asin((hfe_min+hc(end)/2)/(r-pontT(end)+hc(end)/2));
        xcbar = (r-pontT-hc/2).*cos(alpha_pont);
        ycbar = (r-pontT-hc/2).*sin(alpha_pont);
        xpont = (r-pontT).*cos(alpha);
        ypont = (r-pontT).*sin(alpha);
    end
end

% update hc in geo
geo.hc = hc;
geo.hc_pu = hc_pu;

%% barrier angle check
Bx0 = xcbar-(ycbar./tan(beta));
B1k = Bx0-hc/2./sin(beta);
B2k = Bx0+hc/2./sin(beta);

%%MODIFICA

% intShaft=B1k(end)-Ar-hfe_min;
% if intShaft<0
%     disp('#2 barrier cross shaft, angle limit')
%     tmp = calc_Vtype_angle_lim(geo);
%     beta(end) = tmp(end);
%     Bx0   = xcbar-(ycbar./tan(beta));
%     B1k   = Bx0-hc/2./sin(beta);
%     B2k   = Bx0+hc/2./sin(beta);
% end

beta_pu = (1-beta*2/pi)*p/(p-1);

geo.betaPMshape = beta_pu;
geo.beta    = beta;

betaLim = calc_Vtype_angle_lim(geo);
betaLim(1:end-1) = pi/2/p;

%% Check barrier intersections

for ii=2:nlay
    if B2k(ii)>(B1k(ii-1)-hfe_min)
        disp('#2 barriers intersection')
        B2k(ii) = B1k(ii-1)-hfe_min;
        
        aTmp = (B2k(ii)-xcbar(ii)).^2-(hc(ii)/2).^2;
        bTmp = -2*(B2k(ii)-xcbar(ii)).*(0-ycbar(ii));
        cTmp = (0-ycbar(ii)).^2-(hc(ii)/2).^2;
        m1 = (-bTmp+(bTmp.^2-4*aTmp.*cTmp).^0.5)./(2*aTmp);
        m2 = (-bTmp-(bTmp.^2-4*aTmp.*cTmp).^0.5)./(2*aTmp);
        
        beta1 = atan(m1);
        beta2 = atan(m2);
        
        if beta1<0
            beta1=beta1+pi;
        end
        if beta1>pi/2
            beta1 = pi/2;
        end
%         if beta1<betaLim(ii)                  %MODIFICA               
%             beta1 = betaLim(ii);
%         end
        
        if beta2<0
            beta2=beta2+pi;
        end
        if beta2>pi/2
            beta2 = pi/2;
        end
%         if beta2<betaLim(ii)                  %MODIFICA           
%             beta2 = betaLim(ii);
%         end
        
        if beta1<beta2
            beta(ii)=beta1;
        else
            beta(ii)=beta2;
        end
        
        Bx0(ii) = xcbar(ii)-(ycbar(ii)./tan(beta(ii)));
        B1k(ii) = Bx0(ii)-hc(ii)/2./sin(beta(ii));
        B2k(ii) = Bx0(ii)+hc(ii)/2./sin(beta(ii));
        
    end
end

beta_pu = (1-beta*2/pi)*p/(p-1);

%% barrier intersection (more advanced)

for ii=2:nlay
    
    [a2,b2,c2] = retta_per_2pti(Bx0(ii),0,xcbar(ii),ycbar(ii));
    m2=-a2/b2;
    q2=-c2/b2;
    [xTmp,yTmp] = intersezione_retta_circonferenza(xcbar(ii-1),ycbar(ii-1),hc(ii-1)/2+hfe_min+hc(ii)/2,m2,q2);
    if imag(xTmp*yTmp)==0
         disp('#3 barriers intersection')
        aTmp = 2*xcbar(ii)*xcbar(ii-1)-xcbar(ii)^2-xcbar(ii-1)^2+(hc(ii-1)/2+hfe_min+hc(ii)/2)^2;
        bTmp = 2*xcbar(ii)*(ycbar(ii)-ycbar(ii-1))-2*xcbar(ii-1)*(ycbar(ii)-ycbar(ii-1));
        cTmp = -(ycbar(ii)-ycbar(ii-1))^2+(hc(ii-1)/2+hfe_min+hc(ii)/2)^2;
        
        m1 = (-bTmp+(bTmp.^2-4*aTmp.*cTmp).^0.5)./(2*aTmp);
        m2 = (-bTmp-(bTmp.^2-4*aTmp.*cTmp).^0.5)./(2*aTmp);
        
        beta1 = atan(m1);
        beta2 = atan(m2);
        
        if beta1<0
            beta1=beta1+pi;
        end
        if beta1>pi/2
            beta1 = pi/2;
        end
%         if beta1<betaLim(ii)           % MODIFICA              
%             beta1 = betaLim(ii);
%         end
        
        if beta2<0
            beta2=beta2+pi;
        end
        if beta2>pi/2
            beta2 = pi/2;
        end
%         if beta2<betaLim(ii)             % MODIFICA   
%             beta2 = betaLim(ii);
%         end
        
        if beta1<beta2
            beta(ii)=beta1;
        else
            beta(ii)=beta2;
        end
        
        Bx0(ii) = xcbar(ii)-(ycbar(ii)./tan(beta(ii)));
        B1k(ii) = Bx0(ii)-hc(ii)/2./sin(beta(ii));
        B2k(ii) = Bx0(ii)+hc(ii)/2./sin(beta(ii));
    end
end

beta_pu = (1-beta*2/pi)*p/(p-1);

%%MODIFICA

%% MODIFICA
%Shaft check #3
%First step
%Check beta_pu #2° version

drtmp=0;
for ii=2:nlay
    if beta_pu(ii)>1
        disp('#4 beta limit')
        beta_pu(ii) = 1;
        beta(ii) = pi/2/p;
        m = tan(beta(ii));
        q = ycbar(ii)-m*xcbar(ii);
        
        Bx0(ii) = -q/m;
        B1k(ii) = Bx0(ii)-hc(ii)/2./sin(beta(ii));
        B2k(ii) = Bx0(ii)+hc(ii)/2./sin(beta(ii));
        
        [xTmp,yTmp] = intersezione_retta_circonferenza(xcbar(ii-1),ycbar(ii-1),hc(ii-1)/2+hfe_min+hc(ii)/2,m,q);
          if imag(xTmp*yTmp)==0
          drtmp = sqrt(xcbar(ii-1)^2+q^2+ycbar(ii-1)^2-2*ycbar(ii-1)*q-(2*m*q-2*xcbar(ii-1)-2*ycbar(ii-1)*m)^2/4/(1+m^2));
          
          hc(ii-1) = 2*(drtmp-hc(ii)/2-hfe_min);
          
          B1k(ii-1) = Bx0(ii-1)-hc(ii-1)/2./sin(beta(ii-1));
          B2k(ii-1) = Bx0(ii)+hc(ii-1)/2./sin(beta(ii-1));
          end
    end
end

hc_pu = hc/g;
geo.hc = hc;
geo.hc_pu = hc_pu;

%% Check for minimum beta, for FEMM mesh compatibility
betaFEMM = (1-(89.5*pi/180)*2/pi)*p/(p-1);
beta_pu(beta_pu<betaFEMM) = 0;
beta = pi/2-beta_pu*(pi/2-pi/2/p);

geo.beta_pu = beta_pu;
geo.beta    = beta;

Bx0 = xcbar-(ycbar./tan(beta));       
B1k = Bx0-hc/2./sin(beta);
B2k = Bx0+hc/2./sin(beta);

%% End barrier points computation
% Circular (for PM max size) - y=m1*x+q1: barrier axes, y=m2*x+q2: normal 
m1 = (ycbar)./(xcbar-Bx0);
q1 = -m1.*Bx0;
m2 = -1./m1;
q2 = ycbar-m2.*xcbar;
xxD1k = xcbar-abs(hc/2.*cos(atan(m2)));
yyD1k = xxD1k.*m2+q2;
xxD2k = xcbar+abs(hc/2.*cos(atan(m2)));
yyD2k = xxD2k.*m2+q2;

% Computation for the constant thickness ribs
ang1 = atan2(yyD1k,xxD1k);
ang2 = atan2(yyD2k,xxD2k);

xTraf1 = (r-pontT).*cos(ang1);
yTraf1 = (r-pontT).*sin(ang1);
xTraf2 = (r-pontT).*cos(ang2);
yTraf2 = (r-pontT).*sin(ang2);

% if |dx|<0.5, the end-barrier is circular and these points are useless
xTraf1(abs(dx)<0.5) = NaN;
yTraf1(abs(dx)<0.5) = NaN;
xTraf2(abs(dx)<0.5) = NaN;
yTraf2(abs(dx)<0.5) = NaN;

%% Radial ribs computation
% Approximation and hypotesis:
% - neglect the effect of the tangential ribs
% - neglect the barrier presence: PM and iron have almost the same density

pontTmp = zeros(1,nlay);
for ii=1:nlay
    X = [Bx0(ii) xpont(ii) r*cos(linspace(alpha_pont(ii),0,21)) Bx0(ii)];
    Y = [0       ypont(ii) r*sin(linspace(alpha_pont(ii),0,21)) 0      ];
    Afe = polyarea(X,Y);
    rG = centroid(X',Y');
    Mfe = 2*Afe*l*1e-9*rhoFE;
    F = Mfe*rG(1)/1000*(nmax*pi/30)^2;
    pontTmp(ii) = F/(sigmaFe*l); % mm
end

if ~radial_ribs_eval
    pontR = pontTmp;
end

pontR(pontR<pont0) = 0;

XpontRadBarSx = zeros(1,nlay);
YpontRadBarSx = zeros(1,nlay);
XpontRadBarDx = zeros(1,nlay);
YpontRadBarDx = zeros(1,nlay);
XpontRadDx    = zeros(1,nlay);
YpontRadDx    = zeros(1,nlay);
XpontRadSx    = zeros(1,nlay);
YpontRadSx    = zeros(1,nlay);
xS01          = zeros(1,nlay);
yS01          = zeros(1,nlay);
xS02          = zeros(1,nlay);
yS02          = zeros(1,nlay);




%check rotor fillet rad
RotorFillet2(RotorFillet2>hc/2)=round(hc((RotorFillet2>hc/2))/2);
RotorFillet1(RotorFillet1<pont0)=pont0;
RotorFillet1(RotorFillet1>hc/2)=round(hc((RotorFillet1>hc/2))/2);
RotorFillet1(RotorFillet1<pont0)=pont0;

for ii=1:nlay % controlla se utilizzare calcolo matriciale invece di ciclo for
    if pontR(ii)==0
        XpontRadDx(ii) = NaN;
        YpontRadDx(ii) = 0;
        XpontRadSx(ii) = NaN;
        YpontRadSx(ii) = 0;
        XpontRadBarDx(ii) = B2k(ii);
        YpontRadBarDx(ii) = 0;
        XpontRadBarSx(ii) = B1k(ii);
        YpontRadBarSx(ii) = 0;
    else
        a0 = 0;
        b0 = 1;
        c0 = -pontR(ii)/2;
        [a1,b1,c1] = retta_per_2pti(B1k(ii),0,xxD1k(ii),yyD1k(ii));
        [a2,b2,c2] = retta_per_2pti(B2k(ii),0,xxD2k(ii),yyD2k(ii));

        [m0,q0,~]=retta_abc2mq(a0,b0,c0);
        [m1,q1,~]=retta_abc2mq(a1,b1,c1);
        [m2,q2,~]=retta_abc2mq(a2,b2,c2);

        yp1=pontR(ii)/2;
        xp1=(yp1-q1)./m1;

        yp2=yp1;
        xp2=(yp2-q2)./m2;

        XpontRadSx(ii)=xp1+RotorFillet1(ii);
        YpontRadSx(ii)=yp1;

        XpontRadDx(ii)=xp2-RotorFillet2(ii);
        YpontRadDx(ii)=yp2;

        a=(1+m1.*m1);
        b=-2*xp1+2*m1.*(q1-yp1);
        c=-RotorFillet1(ii).*RotorFillet1(ii)+xp1.^2+(q1-yp1).^2;
        XpontRadBarSx(ii)=(-b+sqrt(b.^2-4*a.*c))./(2*a);
        YpontRadBarSx(ii)=m1.*XpontRadBarSx(1,(ii))+q1;

        a=(1+m2.*m2);
        b=-2*xp2+2*m2.*(q2-yp2);
        c=-RotorFillet2(ii).*RotorFillet2(ii)+xp2.^2+(q2-yp2).^2;
        XpontRadBarDx(ii)=(-b+sqrt(b.^2-4*a.*c))./(2*a);
        YpontRadBarDx(ii)=m2.*XpontRadBarDx(1,(ii))+q2;

        m1perp=-1/m1;
        q1perp=YpontRadBarSx(ii)-m1perp*XpontRadBarSx(ii);

        m2perp=-1/m2;
        q2perp=YpontRadBarDx(ii)-m2perp*XpontRadBarDx(ii);

        xS01(ii)=XpontRadSx(ii);
        yS01(ii)=m1perp.*xS01(ii)+q1perp;

        xS02(ii)=XpontRadDx(ii);
        yS02(ii)=m2perp.*xS02(ii)+q2perp;
    end
end

X1 = XpontRadBarDx-hc.*sin(beta);
Y1 = YpontRadBarDx+hc.*cos(beta);
X2 = XpontRadBarDx;
Y2 = YpontRadBarDx;

%% MODIFICA 
%Shaft check #3
%Second step
%Cut through barriers (Without pontRad)

FlagCutMag = zeros(1:nlay);
for ii=1:nlay
    [a,b,c] = retta_per_2pti(B1k(ii),0,xxD1k(ii),yyD1k(ii));
    m = -a/b;
    q = -c/b;
    [xTmp,yTmp] = intersezione_retta_circonferenza(0,0,Ar+hfe_min,m,q);
        if (imag(xTmp*yTmp)==0 && xTmp*yTmp >0)
            disp('#5 barriers cut')
            FlagCutMag(ii) = ii;
            XpontRadBarSx(ii) = xTmp+hc(ii)/2*cos(beta(ii));
            YpontRadBarSx(ii) = yTmp+hc(ii)/2*sin(beta(ii));
            XpontRadBarDx(ii) = XpontRadBarSx(ii)+hc(ii)*sin(beta(ii));
            YpontRadBarDx(ii) = YpontRadBarSx(ii)-hc(ii)*cos(beta(ii));
            
            XpontRadSx(ii) = xTmp+pont0*sin(beta(ii));
            YpontRadSx(ii) = yTmp-pont0*cos(beta(ii));
            XpontRadDx(ii) = XpontRadSx(ii)+(hc(ii)-2*pont0)*sin(beta(ii));
            YpontRadDx(ii) = YpontRadSx(ii)-(hc(ii)-2*pont0)*cos(beta(ii));
           
            if XpontRadDx(ii)<XpontRadSx(ii)
                XpontTmp = XpontRadDx(ii);
                YpontTmp = YpontRadDx(ii);
                XpontRadDx(ii) = XpontRadSx(ii);
                YpontRadDx(ii) = YpontRadSx(ii);
                XpontRadSx(ii) = XpontTmp;
                YpontRadSx(ii) = YpontTmp;
            end
            
            if YpontRadBarDx(ii)<pont0
                YpontRadBarDx(ii) = 0;
                m = tan(beta(ii));
                q = -m*xxD2k(ii)+yyD2k(ii);
                XpontRadBarDx(ii) = -q/m;

                YpontRadDx(ii) = YpontRadBarDx(ii);
                XpontRadDx(ii) = XpontRadBarDx(ii);
            end
        
            if YpontRadDx(ii)<pont0
                YpontRadDx(ii) = 0;
                m = -1/tan(beta(ii));
                q = -m*XpontRadSx(ii)+YpontRadSx(ii);
                XpontRadDx(ii) = -q/m;

                YpontRadBarDx(ii) = YpontRadDx(ii);
                m = tan(beta(ii));
                q = -m*xxD2k(ii)+yyD2k(ii);
                XpontRadBarDx(ii) = -q/m;
         
            	if XpontRadBarDx(ii)<XpontRadDx(ii)
                    XpontRadDx(ii) = XpontRadBarDx(ii);
            	end
            end
            
            if YpontRadSx(ii)<pont0
                YpontRadSx(ii) = YpontRadBarSx(ii);
                XpontRadSx(ii) = XpontRadBarSx(ii);
            end       
        end
end

if sum(find(FlagCutMag))>0
X1 = XpontRadBarDx-hc.*sin(beta);
Y1 = YpontRadBarDx+hc.*cos(beta);
X2 = XpontRadBarDx;
Y2 = YpontRadBarDx;
    for ii=1:nlay
        if X1(ii)<XpontRadBarSx(ii)
            X1(ii) = XpontRadBarSx(ii);
            Y1(ii) = YpontRadBarSx(ii);
           	X2(ii) = X1(ii)+hc(ii)*sin(beta(ii));
            Y2(ii) = Y1(ii)-hc(ii)*cos(beta(ii));
        end
    end
end

%Cut through barriers (With pontRad)??
%In the range of speed evaluated, if the barrieres are cut there will be
%not pont Rad.?

%% PM dimensions

b1 = ((X1-xxD1k).^2+(Y1-yyD1k).^2).^0.5;
b2 = ((X2-xxD2k).^2+(Y2-yyD2k).^2).^0.5;

b = zeros(1,nlay);
b(b1<b2) = b1(b1<b2);
b(b2<=b1) = b2(b2<=b1);

AreaCMax = b.*hc;
if min(PMdimC)<0 % per-unit input (during optimization)
    PMdimC = -PMdimC.*b;
end

PMdimC(PMdimC>b) = b(PMdimC>b);
PMdimC = floor(PMdimC*100)/100;
PMdimC(PMdimC<0) = 0;

AreaC = PMdimC.*hc;

% Check clearance size: the minimum PM thickness is set to pont0
PMclearLim = hc-pont0*ones(size(hc));
PMclearC(PMclearC>PMclearLim) = PMclearLim(PMclearC>PMclearLim);

xPMC1b = X1-PMclearC.*cos(beta+pi/2);
yPMC1b = Y1-PMclearC.*sin(beta+pi/2);
xPMC2b = X2;
yPMC2b = Y2;

xPMC1t = xPMC1b+PMdimC.*cos(beta);
yPMC1t = yPMC1b+PMdimC.*sin(beta);
xPMC2t = xPMC2b+PMdimC.*cos(beta);
yPMC2t = yPMC2b+PMdimC.*sin(beta);

% second area (for compatibility with Seg and circular)
AreaE    = zeros(1,nlay);
AreaEMax = zeros(1,nlay);

xPME1b = zeros(1,nlay);
yPME1b = zeros(1,nlay);
xPME2b = zeros(1,nlay);
yPME2b = zeros(1,nlay);
xPME1t = zeros(1,nlay);
yPME1t = zeros(1,nlay);
xPME2t = zeros(1,nlay);
yPME2t = zeros(1,nlay);

PMdimE = zeros(size(PMdimC));

%% Block labels 
xc = (xPMC1b+xPMC1t+xPMC2b+xPMC2t)/4;
yc = (yPMC1b+yPMC1t+yPMC2b+yPMC2t)/4;
xmag = -sin(beta).*ones(1,nlay);
ymag = +cos(beta).*ones(1,nlay);

xair = zeros(1,3*nlay);
yair = zeros(1,3*nlay);
indexZone = repmat([1,2,3],[1,nlay]);
% air zones:
% 1 --> between the PM and the d-axis, in the triangular air region
% 2 --> close to the radial ribs, inside its boundaries (as the other geometries)
% 3 --> between PM and tangential radial rib

%MODIFICA After barriers cut

for ii=0:nlay-1
    if FlagCutMag(ii+1)>0
       xair(1+ii*3) = (xPMC1b(ii+1)+xPMC2b(ii+1)+XpontRadSx(ii+1)+XpontRadDx(ii+1))/4; % air between the PM and the d-axis
       yair(1+ii*3) = (yPMC1b(ii+1)+yPMC2b(ii+1)+YpontRadSx(ii+1)+YpontRadDx(ii+1))/4;
    else
       xair(1+ii*3) = (xPMC1b(ii+1)+xPMC2b(ii+1)+XpontRadBarSx(ii+1)+XpontRadBarDx(ii+1))/4; % air between the PM and the d-axis
       yair(1+ii*3) = (yPMC1b(ii+1)+yPMC2b(ii+1)+YpontRadBarSx(ii+1)+YpontRadBarDx(ii+1))/4; 
    end
end
xair(indexZone==2) = (XpontRadBarDx+XpontRadBarSx+XpontRadDx+XpontRadSx)/4;   % air of the radial rib
yair(indexZone==2) = (YpontRadBarDx+YpontRadBarSx+YpontRadDx+YpontRadSx)/4;
xair(indexZone==3) = (xxD1k+xxD2k+xpont)/3; % air between the PM and the tangential rib
yair(indexZone==3) = (yyD1k+yyD2k+ypont)/3;

% filt xc and xair

% Rules:
% PM zone: deactivated if PM not present
% air zone 3: always activated
% air zone 1: activated only if beta>0 and PM present
% air sone 2: activated only if beta==0 and PM present

xc(PMdimC==0) = NaN;

xair(indexZone==1) = xair(indexZone==1).*(PMdimC./abs(PMdimC));
xair(indexZone==1) = xair(indexZone==1).*(beta_pu./abs(beta_pu));

xair(indexZone==2) = xair(indexZone==2).*(PMdimC./abs(PMdimC));
fTmp = nan(1,nlay);
fTmp(beta_pu==0) = 1;
xair(indexZone==2) = xair(indexZone==2).*fTmp;
filtClear = ones(size(PMclearC));
filtClear(PMclearC>0) = NaN;
xair(indexZone==1) = xair(indexZone==1).*filtClear; % if PMC are with clearance, rad rib labels must be disabled

% MODIFICA End Barrier information for block label 

for ii=nlay:-1:1
    if (Ar>B1k(ii))
        B1k(ii) = B1k(ii)+(Ar-B1k(ii))*1.1;
    end
        
    if FlagCutMag(ii)>0
        B1k(end) = B1k(ii);
    end
end

%% Output data
geo.pontR       = pontR;
geo.hc          = hc;
geo.hc_pu       = hc_pu;
geo.betaPMshape = beta_pu;
geo.beta        = beta;
geo.B1k         = B1k;
geo.B2k         = B2k;
geo.Bx0         = Bx0;
geo.xxD1k       = xxD1k;
geo.yyD1k       = yyD1k;
geo.xxD2k       = xxD2k;
geo.yyD2k       = yyD2k;
geo.xPMC1b      = xPMC1b;
geo.yPMC1b      = yPMC1b;
geo.yPMC2b      = yPMC2b;
% geo.xPMC2t      = xPMC2t;
% geo.yPMC2t      = yPMC2t;

geo.RotorFillet1 = RotorFillet1;
geo.RotorFillet2 = RotorFillet2;


geo.AreaCMax  = AreaCMax;
geo.AreaEMax  = AreaEMax;
geo.AreaC     = AreaC;
geo.AreaE     = AreaE;
geo.PMdim     = [PMdimC;PMdimE];
geo.PMclear   = [PMclearC;ones(size(PMclearC))];

% end barrier
temp.B1k        = B1k;
temp.B2k        = B2k;
temp.Bx0        = Bx0;
temp.xpont      = xpont;
temp.ypont      = ypont;
temp.xxD1k      = xxD1k;
temp.yyD1k      = yyD1k;
temp.xxD2k      = xxD2k;
temp.yyD2k      = yyD2k;
temp.XcRibTraf1 = xcbar;
temp.YcRibTraf1 = ycbar;
temp.XcRibTraf2 = xcbar;
temp.YcRibTraf2 = ycbar;
temp.RcRibTraf2 = hc/2;
temp.xTraf1     = xTraf1;
temp.yTraf1     = yTraf1;
temp.xTraf2     = xTraf2;
temp.yTraf2     = yTraf2;

% labels
temp.xc       = xc;
temp.yc       = yc;
temp.xair     = xair;
temp.yair     = yair;
temp.xmag     = xmag;
temp.ymag     = ymag;
temp.zmag     = zeros(size(xmag));
 
% radial ribs
temp.XpontRadDx      = XpontRadDx;
temp.YpontRadDx      = YpontRadDx;
temp.XpontRadSx      = XpontRadSx;
temp.YpontRadSx      = YpontRadSx;
temp.XpontRadBarDx   = XpontRadBarDx;
temp.XpontRadBarSx   = XpontRadBarSx;
temp.YpontRadBarDx   = YpontRadBarDx;
temp.YpontRadBarSx   = YpontRadBarSx;
temp.xS01            = xS01;
temp.yS01            = yS01;
temp.xS02            = xS02;
temp.yS02            = yS02;
% temp.XcRaccpontRadSx = XcRaccpontRadSx;
% temp.YcRaccpontRadSx = YcRaccpontRadSx;
% temp.RcRaccpontRadSx = RcRaccpontRadSx;
% temp.XcRaccpontRadDx = XcRaccpontRadDx;
% temp.YcRaccpontRadDx = YcRaccpontRadDx;
% temp.RcRaccpontRadDx = RcRaccpontRadDx;

% magnet
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

temp.mirrorFlag    = ones(size(xc));
temp.mirrorFlagAir = ones(size(xair));

%mat.LayerMag.Br = [mat.LayerMag.Br mat.LayerMag.Br];
mat.LayerMag.Br = mat.LayerMag.Br(1).*ones(1,nlay*2);

