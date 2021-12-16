% Copyright 2014
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

function [geo,mat,temp]=nodes_rotor_Circ_dx(geo,mat)

r = geo.r;                      % Raggio del rotore al traferro
x0 = geo.x0;                    % Centro fittizio
rshaft = geo.Ar;                % Raggio albero
Ar=geo.Ar;
l = geo.l;                      % Lunghezza pacco
g = geo.g;                      % Traferro
pont0 = geo.pont0;              % minimum mechanical tolerance
pontT = geo.pontT;              % Airgap ribs [mm]

p = geo.p;                      % Paia poli
nlay = geo.nlay;                % N° layers

dalpha = geo.dalpha;            % Angoli dalpha
% Eval alpha
alpha = cumsum(dalpha);
dx=geo.dx;

racc_pont = geo.racc_pont;      % racc_pont=1*pont0 <- per i ponticelli radiali.
% ang_pont0 = geo.ang_pont0;      % Ampiezza dell'angolo (in gradi) da spazzare con  raggio r in modo da ottenre un arco lungo pont0

nmax = geo.nmax;                % Velocità max (rpm) per la valutazione della sollecitazione centrifuga più gravosa (-> ponticelli)
hfe_min=geo.hfe_min;
sigma_max = mat.Rotor.sigma_max;    % snervamento materiale [MPa]
rhoFE = mat.Rotor.kgm3;             % densità del ferro di rotore [kg/m3]
rhoPM = mat.LayerMag.kgm3;          % densità magneti [kg/m3]

RotorFilletTan1 = geo.RotorFilletTan1;
RotorFilletTan2 = geo.RotorFilletTan2;
RotorFilletTan1(~isfinite(RotorFilletTan2)) =  NaN;
RotorFilletTan2(~isfinite(RotorFilletTan1)) =  NaN;


XcRibTraf1=zeros(1,nlay);
XcRibTraf2=zeros(1,nlay);
YcRibTraf1=zeros(1,nlay);
YcRibTraf2=zeros(1,nlay);
xxD1k=zeros(1,nlay);
yyD1k=zeros(1,nlay);
xxD2k=zeros(1,nlay);
yyD2k=zeros(1,nlay);



% CENTRO FITTIZIO DI COORDINATE (x0,0)
beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,r,x0);
% La funzione calc_apertura_cerchio riceve in input le coordinate polari,
% (r, alpha) (alpha in rad), di un punto generico. Queste sono calcolate
% rispetto al centro (0,0). In output la funzione restituisce l'apertura
% angolare (in rad) dello stesso punto rispetto al centro preso come
% riferimento (ha coordinate:(x0,0)).
% I punti di cui, in questo caso, si calcolano le aperture angolari rispetto
% al centro di riferimento sono i punti mediani delle barriere, presi in
% corrispondenza del traferro.

rbeta = (x0 - r * cos(alpha*pi/180))./(cos(beta*pi/180)); 
% Di questi stessi punti, si calcolano anche le distanze dal centro (x0,0)
% e le si memorizzano nel vettore rbeta.

[xpont,ypont] = calc_intersezione_cerchi(r-pontT, rbeta, x0);

LowDimBarrier=zeros(1,nlay);
for ii=1:nlay
    if (not(isreal(xpont(ii)))||not(isreal(ypont(ii))))
        xpont(ii)=(r-2*pontT(ii)).*cos(alpha(ii)*pi/180);
        ypont(ii)=(r-2*pontT(ii)).*sin(alpha(ii)*pi/180);
        LowDimBarrier(ii)=1;
    end
end

rpont_x0=sqrt(ypont.^2+(x0-xpont).^2);
[alphapont,rpont] = cart2pol(xpont,ypont);
Bx0=x0-(rpont_x0);

% Determination of air thickness and check the feasibility of the geometry
geo.Bx0=Bx0; % Initialization of central non-moved line of the flux barrier
if geo.delta_FBS==0
    geo = calcHcCheckGeoControlwDx(geo);
    B1k=geo.B1k;
    B2k=geo.B2k;
else
    B1k=Bx0-geo.hc/2+geo.dx.*geo.hc/2;
    B2k=Bx0+geo.hc/2+geo.dx.*geo.hc/2;
    geo.B1k=B1k;
    geo.B2k=B2k;
end

hc=B2k-B1k;

ptmp=find(abs(hc)<pont0/2);
B2k(ptmp)=min([B1k(ptmp),B2k(ptmp)])+pont0/2;

% Intersezione circonferenze punti al traferro:
[xTraf2,yTraf2] = calc_intersezione_cerchi(r-pontT, x0-B2k, x0);
[xTraf1,yTraf1] = calc_intersezione_cerchi(r-pontT, x0-B1k, x0);

for ii=1:nlay
    if (not(isreal(xTraf1(ii)))||not(isreal(yTraf1(ii)))||not(isreal(xTraf2(ii)))||not(isreal(yTraf2(ii)))||hc(ii)<0)
        xpont(ii)=(r-2*pontT(ii)).*cos(alpha(ii)*pi/180);
        ypont(ii)=(r-2*pontT(ii)).*sin(alpha(ii)*pi/180);
        LowDimBarrier(ii)=1;
    end
end

% Barriers tips nodes (points 1 and 2 + centers)
RotorFilletTan1 (RotorFilletTan1<pont0*ones(1,nlay) & isfinite(RotorFilletTan1)) = pont0;
RotorFilletTan2 (RotorFilletTan2<pont0*ones(1,nlay) & isfinite(RotorFilletTan2)) = pont0;
RotorFilletTan1 (RotorFilletTan1>geo.hc/2 & isfinite(RotorFilletTan1)) = geo.hc(RotorFilletTan1>geo.hc/2 & isfinite(RotorFilletTan1))/2;
RotorFilletTan2 (RotorFilletTan2>geo.hc/2 & isfinite(RotorFilletTan2)) = geo.hc(RotorFilletTan2>geo.hc/2 & isfinite(RotorFilletTan2))/2;

for ii=1:nlay
    if LowDimBarrier(ii)==1
        xxD1k(ii)=B1k(ii);
        yyD1k(ii)=0;
        xxD2k(ii)=B2k(ii);
        yyD2k(ii)=0;
    else

        if ~(isfinite(RotorFilletTan1(ii)) & isfinite(RotorFilletTan2(ii)))
            [xt,yt,xc,yc,rc]=cir_tg_2cir(xpont(ii),ypont(ii),r-pontT(ii),x0,0,x0-B1k(ii));
            xxD1k(ii)=xt;      % D1: end of outer semi arc
            yyD1k(ii)=yt;
            XcRibTraf1(ii)=xc; % center of semi arc 1
            YcRibTraf1(ii)=yc;
            [xt,yt,xc,yc,rc]=cir_tg_2cir(xpont(ii),ypont(ii),r-pontT(ii),x0,0,x0-B2k(ii));
            xxD2k(ii)=xt;      % D2: end of inner semi arc
            yyD2k(ii)=yt;
            XcRibTraf2(ii)=xc; % center of semi arc 2
            YcRibTraf2(ii)=yc;
            
        else
            [xtemp,ytemp] = calc_intersezione_cerchi(r-pontT(ii)-RotorFilletTan1(ii),x0-B1k(ii)-RotorFilletTan1(ii),x0);
            xcRac1(ii) = real(xtemp);
            ycRac1(ii) = real(ytemp);
            [xtemp,ytemp]=intersezione_cerchi(xcRac1(ii),ycRac1(ii),RotorFilletTan1(ii),x0,0,x0-B1k(ii));
            xxD1k(ii)=real(xtemp(2));
            yyD1k(ii)=real(ytemp(2));
            [xtemp,ytemp]=intersezione_cerchi(xcRac1(ii),ycRac1(ii),RotorFilletTan1(ii),0,0,r-pontT(ii));
            xxE1k(ii)=real(xtemp(2));
            yyE1k(ii)=real(ytemp(2));

            [xtemp,ytemp] = calc_intersezione_cerchi(r-pontT(ii)-RotorFilletTan2(ii),x0-B2k(ii)+RotorFilletTan2(ii),x0);
            xcRac2(ii) = real(xtemp);
            ycRac2(ii) = real(ytemp);
            [xtemp,ytemp]=intersezione_cerchi(xcRac2(ii),ycRac2(ii),RotorFilletTan2(ii),x0,0,x0-B2k(ii));
            xxD2k(ii)=real(xtemp(2));
            yyD2k(ii)=real(ytemp(2));
            [xtemp,ytemp]=intersezione_cerchi(xcRac2(ii),ycRac2(ii),RotorFilletTan2(ii),0,0,r-pontT(ii));
            xxE2k(ii)=real(xtemp(2));
            yyE2k(ii)=real(ytemp(2));

            temp.xcRac1 = xcRac1;
            temp.ycRac1 = ycRac1;
            temp.xxE1k = xxE1k;
            temp.yyE1k = yyE1k;
            temp.xcRac2 = xcRac2;
            temp.ycRac2 = ycRac2;
            temp.xxE2k = xxE2k;
            temp.yyE2k = yyE2k;

        end
    end
end

temp.B1k=B1k;
temp.B2k=B2k;
temp.Bx0=Bx0;

[temp,geo] = calc_ribs_rad_fun(geo,mat,temp);

% Points for radial ribs
XpontRadDx    = temp.XpontRadDx;
YpontRadDx    = temp.YpontRadDx;
XpontRadSx    = temp.XpontRadSx;
YpontRadSx    = temp.YpontRadSx;
XpontRadBarDx = temp.XpontRadBarDx;
XpontRadBarSx = temp.XpontRadBarSx;
YpontRadBarDx = temp.YpontRadBarDx;
YpontRadBarSx = temp.YpontRadBarSx;

temp.xc=(xxD1k+xxD2k)/2;
temp.yc=(yyD1k+yyD2k)/2;

% Determining  Magnet Area
temp.xpont=xpont;
temp.ypont=ypont;
temp.xxD1k=xxD1k;
temp.yyD1k=yyD1k;
temp.xxD2k=xxD2k;
temp.yyD2k=yyD2k;
temp.XcRibTraf1=XcRibTraf1;
temp.YcRibTraf1=YcRibTraf1;
temp.XcRibTraf2=XcRibTraf2;
temp.YcRibTraf2=YcRibTraf2;

[geo,mat,temp] = PMdefinition_Circ(geo,mat,temp);

% add control if the construction of barrier is too small, so barrier is
% constructed like diamond
temp.LowDimBarrier=LowDimBarrier;
hf=[r,B1k]-[B2k,Ar]; %calcolo dei Delta fi ferro di rotore
geo.hf = hf;
% geo.pont = pont;
geo.hc=hc;
geo.xpont=xpont;
geo.ypont=ypont;

% barrier transverse dimension (for permeance evaluation)
temp_r_arcs = (geo.r_all(1:2:end-1)+geo.r_all(2:2:end))/2;
temp_x_arc_ends = x0 - (xpont);
% equivalent sk = arc barrier, disregarding the end barrier radius
sk = temp_r_arcs.*acos(temp_x_arc_ends./temp_r_arcs);

geo.sk = sk;
geo.pbk = geo.sk ./ geo.hc;
geo.la = sum(geo.hc)/geo.r;
geo.lfe = sum(geo.hf)/geo.r;
geo.ly = (geo.R - (geo.r + geo.g + geo.lt))/geo.r;
geo.B1k=B1k;
geo.B2k=B2k;
geo.xxD1k=xxD1k;
geo.yyD1k=yyD1k;
geo.xxD2k=xxD2k;
geo.yyD2k=yyD2k;
geo.hc=hc;
geo.RotorFilletTan1 = RotorFilletTan1;
geo.RotorFilletTan2 = RotorFilletTan2;

