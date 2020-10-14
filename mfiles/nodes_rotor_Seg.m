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

r = geo.r;                      % Raggio del rotore al traferro
x0 = geo.x0;                    % Centro fittizio
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

delta_FBS = geo.delta_FBS;      % flux barrier shift angle

racc_pont = geo.racc_pont;      % racc_pont=1*pont0 <- per i ponticelli radiali.
ang_pont0 = geo.ang_pont0;      % Ampiezza dell'angolo (in gradi) da spazzare con  raggio r in modo da ottenre un arco lungo pont0

nmax = geo.nmax;                % Velocità max (rpm) per la valutazione della sollecitazione centrifuga più gravosa (-> ponticelli)
hfe_min=geo.hfe_min;
sigma_max = mat.Rotor.sigma_max;    % snervamento materiale [MPa]
rhoFE = mat.Rotor.kgm3;             % densità del ferro di rotore [kg/m3]
rhoPM = mat.LayerMag.kgm3;          % densità magneti [kg/m3]

alphaRot=alpha*pi/180+pi/2/p;   % angoli barriere di flux al traferro meccanici

rCirRib=1*pont0;
rlim=r;

RotorFillet=geo.RotorFillet;
% DTrasl=geo.DTrasl;

arcLayTraf1=zeros(1,nlay);
arcLayTraf2=zeros(1,nlay);
XcRibTraf1=zeros(1,nlay);
XcRibTraf2=zeros(1,nlay);
YcRibTraf1=zeros(1,nlay);
YcRibTraf2=zeros(1,nlay);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ampiezza di un arco lungo pont0
% ang_pont0 = pont0 / r * 180/pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISEGNO DELLE BARRIERE E DEI PONTICELLI
% INTRODUZIONE DI ALCUNE GRANDEZZE GEOMETRICHE UTILI E RIFERIMENTO AD UN CENTRO FITTIZIO DI COORDINATE (x0,0)

beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,r,x0);                  % La funzione calc_apertura_cerchio riceve in input le coordinate polari,
% (r, alpha) (alpha in rad), di un punto generico. Queste sono calcolate
% rispetto al centro (0,0). In output la funzione restituisce l'apertura
% angolare (in rad) dello stesso punto rispetto al centro preso come
% riferimento (ha coordinate:(x0,0)).
% I punti di cui, in questo caso, si calcolano le aperture angolari rispetto
% al centro di riferimento sono i punti mediani delle barriere, presi in
% corrispondenza del traferro.

rbeta = (x0 - r * cos(alpha*pi/180))./(cos(beta*pi/180));                      % Di questi stessi punti, si calcolano anche le distanze dal centro (x0,0)
% e le si memorizzano nel vettore rbeta.

[xpont,ypont] = calc_intersezione_cerchi(r-pontT, rbeta, x0);

LowDimBarrier=zeros(1,nlay);
for ii=1:nlay
    if (not(isreal(xpont(ii)))||not(isreal(ypont(ii))))
        xpont(ii)=(r-2*pontT(ii))*cos(alpha(ii)*pi/180);
        ypont(ii)=(r-2*pontT(ii))*sin(alpha(ii)*pi/180);
        LowDimBarrier(ii)=1;
    end
end

rpont_x0=sqrt(ypont.^2+(x0-xpont).^2);
[alphapont,rpont] = cart2pol(xpont,ypont);
Bx0=x0-(rpont_x0); geo.Bx0=Bx0;

% Determination of air thickness and check the feasibility of the geometry
geo.Bx0=Bx0; % Initialization of central non-moved line of the flux barrier
if delta_FBS==0
    geo = calcHcCheckGeoControlwDx(geo);
    B1k=geo.B1k;
    B2k=geo.B2k;
else
    B1k=Bx0-geo.hc/2+geo.dx.*geo.hc/2;
    B2k=Bx0+geo.hc/2+geo.dx.*geo.hc/2;
    geo.B1k=B1k;
    geo.B2k=B2k;
end
% 
% B1k(3) = 25.6540;
% B2k(3) = 33.5280;

hc=B2k-B1k;

ptmp=find(abs(hc)<pont0/2);
B2k(ptmp)=min([B1k(ptmp),B2k(ptmp)])+pont0/2;
clear ptmp;

% Intersezione circonferenze punti al traferro:
% [xTraf2,yTraf2] = calc_intersezione_cerchi(r-pontT, x0-B2k, x0);
% [xTraf1,yTraf1] = calc_intersezione_cerchi(r-pontT, x0-B1k, x0);

m0 = tan(pi/2/p+delta_FBS/2)*ones(1,nlay);
q0 = ypont-m0.*xpont;

x2 = xpont+(B2k-Bx0).*sin(pi/2/p+delta_FBS/2);
y2 = ypont-(B2k-Bx0).*cos(pi/2/p+delta_FBS/2);
m2 = m0;
q2 = y2-m2.*x2;
[xTraf2,yTraf2] = intersezione_retta_circonferenza(0,0,(r-pontT),m2,q2);

x1 = xpont+(B1k-Bx0).*sin(pi/2/p+delta_FBS/2);
y1 = ypont-(B1k-Bx0).*cos(pi/2/p+delta_FBS/2);
m1 = m0;
q1 = y1-m1.*x1;
[xTraf1,yTraf1] = intersezione_retta_circonferenza(0,0,(r-pontT),m1,q1);


for ii=1:nlay
    if (not(isreal(xTraf1(ii)))||not(isreal(yTraf1(ii)))||not(isreal(xTraf2(ii)))||not(isreal(yTraf2(ii)))||hc(ii)<0)
        xpont(ii)=(r-2*pontT(ii))*cos(alpha(ii)*pi/180);
        ypont(ii)=(r-2*pontT(ii))*sin(alpha(ii)*pi/180);
        LowDimBarrier(ii)=1;
    end
end


% 1° barriera di flux: la prima barrira di flusso presenta raccordo
% tondeggiante al traferro diviso in parti uguali pari alla più piccola
% delle distanze da pont0.
% dhoriz=calc_distanza_punti([B2k(1),0],[xTraf2(1),0]);
% dDiag=calc_distanza_punti([B2k(1),0],[xTraf2(1),yTraf2(1)]);
% if (dhoriz<=dDiag/4)
%     xtraf1=xTraf1(1);
%     ytraf1=yTraf1(1);
%     xtraf2=xTraf2(1);
%     ytraf2=yTraf2(1);
%     
%     dt11=abs((xpont(1)-xtraf1)+1j*(ytraf1-ypont(1)));
%     dt12=abs((xpont(1)-xtraf2)+1j*(ytraf2-ypont(1)));
%     [dt1,posMin]=min([dt11,dt12]);
%     m=tan(pi/2/p);
%     B=2*(ypont(1)-m*xpont(1));
%     C=ypont(1)^2+m^2*xpont(1)^2-2*m*xpont(1)*ypont(1)-dt1^2*(1+m^2);
%     q=B/2-sqrt(B^2/4-C);
%     mprimo=-1/m;
%     qprimo=ypont(1)-mprimo*xpont(1);
%     if (posMin==1)
%         xTraf2(1)=(qprimo-q)/(m-mprimo);
%         yTraf2(1)=m*xTraf2(1)+q;
%         xTraf1(1)=xtraf1(1);
%         yTraf1(1)=ytraf1(1);
%     else
%         xTraf2(1)=(qprimo-q)/(m-mprimo);
%         yTraf2(1)=m*xTraf2(1)+q;
%         xTraf1(1)=xtraf1(1);
%         yTraf1(1)=ytraf1(1);
%         
%     end
% end

% Control of the distance between (xTraf1,2;yTraf1,2) to keep minimum iron
for ii=1:nlay-1
    dPto=calc_distanza_punti([xTraf2(ii+1),yTraf2(ii+1)],[xTraf1(ii),yTraf1(ii)]);
    dP1=calc_distanza_punti([xTraf1(ii),yTraf1(ii)],[xpont(ii),ypont(ii)]);
    dP2=calc_distanza_punti([xTraf2(ii+1),yTraf2(ii+1)],[xpont(ii+1),ypont(ii+1)]);
    if(dPto<0.0*hfe_min)
        if (dP1>dP2)
            [a,b,c]=retta_per_2pti(xTraf2(ii+1),yTraf2(ii+1),xTraf1(ii),yTraf1(ii));
            m=-a/b;
            A=1;
            B=-2*xTraf2(ii+1);
            C=xTraf2(ii+1)^2-(1.5*hfe_min)^2/(m^2+1);
            x_temp=roots([A,B,C]);
            xTraf1(ii)=x_temp(1);
            yTraf1(ii)=m*xTraf1(ii)-c/b;
        else
            [a,b,c]=retta_per_2pti(xTraf2(ii+1),yTraf2(ii+1),xTraf1(ii),yTraf1(ii));
            m=-a/b;
            A=1;
            B=-2*xTraf1(ii);
            C=xTraf1(ii)^2-(1.5*hfe_min)^2/(m^2+1);
            x_temp=roots([A,B,C]);
            xTraf2(ii+1)=x_temp(2);
            yTraf2(ii+1)=m*xTraf2(ii+1)-c/b;
        end
    end
end

%
% Intersezione tra rette che compongono i lati delle barriere di flux:
%
% B1k(3) = 25.6540;
% B2k(3) = 33.5280;

% B1k = [55.1180   41.4020   25.6540];
% B2k = [61.7220   47.7520   33.5280];

% retta verticale
for ii=1:nlay
    if LowDimBarrier(ii)==1
        XpBar2(ii)=B2k(ii);
        YpBar2(ii)=0;
    else
        a1=1;
        b1=0;
        c1=-B2k(ii);
        m2=tan(pi/2/p+delta_FBS/2);
        a2=m2;
        b2=-1;
        c2=(yTraf2(ii)-m2*xTraf2(ii));
        [x,y]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
        XpBar2(ii)=x;
        YpBar2(ii)=y;
    end
end

for ii=1:nlay
    if LowDimBarrier(ii)==1
        XpBar1(ii)=B1k(ii);
        YpBar1(ii)=0;
    else
        a1=1;
        b1=0;
        c1=-B1k(ii);
        m2=tan(pi/2/p+delta_FBS/2);
        a2=m2;
        b2=-1;
        c2=(yTraf1(ii)-m2*xTraf1(ii));
        [x,y]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
        XpBar1(ii)=x;
        YpBar1(ii)=y;
    end
end

for ii=1:nlay
    if LowDimBarrier(ii)==1
        xxD1k(ii)=B1k(ii);
        yyD1k(ii)=0;
    else
        [xD1k,yD1k,xc,yc,rc]=tg_cir(XpBar1(ii),YpBar1(ii),xTraf1(ii),yTraf1(ii),xpont(ii),ypont(ii));
        xxD1k(ii)=xD1k;
        yyD1k(ii)=yD1k;
        
        angT2=atan2(yD1k-yc,xD1k-xc);
        angT1=atan2(ypont(ii)-yc,xpont(ii)-xc);
        if (angT2<0)
            angT2=pi+angT2;
        end
        arcLayTraf1(ii)=(angT2-angT1)*180/pi;
        XcRibTraf1(ii)=xc;
        YcRibTraf1(ii)=yc;
        
        clear XCer YCer
        XCer=xc+rc*cos(linspace(angT1,angT2,20));
        YCer=yc+rc*sin(linspace(angT1,angT2,20));
    end
end

for ii=1:nlay
    if LowDimBarrier(ii)==1
        xxD2k(ii)=B2k(ii);
        yyD2k(ii)=0;
    else
        [xD2k,yD2k,xc,yc,rc]=tg_cir(XpBar2(ii),YpBar2(ii),xTraf2(ii),yTraf2(ii),xpont(ii),ypont(ii));
        xxD2k(ii)=xD2k;
        yyD2k(ii)=yD2k;
        
        angT1=atan2(yD2k-yc,xD2k-xc);
        angT2=atan2(ypont(ii)-yc,xpont(ii)-xc);
        XcRibTraf2(ii)=xc;
        YcRibTraf2(ii)=yc;
        arcLayTraf2(ii)=(angT2-angT1)*180/pi;
        xxD2k(ii)=xD2k;
        yyD2k(ii)=yD2k;
        clear XCer YCer
        XCer=xc+rc*cos(linspace(angT1,angT2,20));
        YCer=yc+rc*sin(linspace(angT1,angT2,20));
        %         figure(100);hold on;plot(XCer,YCer,'+m');hold off;
    end
end

% Controllo di sicurezza per evitare che la circonferenza del rib traf termini dopo il lato di barriera
for ii=1:nlay
    if LowDimBarrier(ii)==0
        if (xxD2k(ii)<=XpBar2(ii))
            if(XpBar2(ii)>xpont(ii))
                XpBar2(ii)=xpont(ii);
                B2k(ii)=xpont(ii);
            end
            rc=sqrt((xxD2k(ii)-XcRibTraf2(ii))^2+(yyD2k(ii)-YcRibTraf2(ii))^2);
            
            xxD2k(ii)=XpBar2(ii);
            yyD2k(ii)=-sqrt(rc^2-(XpBar2(ii)-XcRibTraf2(ii))^2)+YcRibTraf2(ii);
            YpBar2(ii)=yyD2k(ii);
            angT1=atan2(yyD2k(ii)-YcRibTraf2(ii),xxD2k(ii)-XcRibTraf2(ii));
            angT2=atan2(ypont(ii)-YcRibTraf2(ii),xpont(ii)-XcRibTraf2(ii));
            arcLayTraf2(ii)=(angT2-angT1)*180/pi;
        end
    end
end
for ii=1:nlay
    if LowDimBarrier(ii)==0
        if (xxD1k(ii)<=XpBar1(ii))
            
            rc=sqrt((xxD1k(ii)-XcRibTraf1(ii))^2+(yyD1k(ii)-YcRibTraf1(ii))^2);
            
            xxD1k(ii)=XpBar1(ii);
            yyD1k(ii)=sqrt(rc^2-(XpBar1(ii)-XcRibTraf1(ii))^2)+YcRibTraf1(ii);
            YpBar1(ii)=yyD1k(ii);
            angT1=atan2(YcRibTraf1(ii)-yyD1k(ii),xxD1k(ii)-XcRibTraf1(ii));
            angT2=atan2(ypont(ii)-YcRibTraf1(ii),xpont(ii)-XcRibTraf1(ii));
            arcLayTraf1(ii)=(angT2-angT1)*180/pi;
        end
    end
end

if(xxD2k(1)<xpont(1) && yyD2k(1)>=ypont(1))
    YpBar2(1)=YpBar2(1)/2;
    yyD2k(1)=YpBar2(1);
end

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raccordo barriere di flux:
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Si esegue un raccordo circolare al posto dello spigolo
XcRacc_B1=NaN*ones(1,nlay);
YcRacc_B1=NaN*ones(1,nlay);
xRaccR1_B1=NaN*ones(1,nlay);
yRaccR1_B1=NaN*ones(1,nlay);
xRaccR2_B1=NaN*ones(1,nlay);
yRaccR2_B1=NaN*ones(1,nlay);

for ii=1:length(B1k)
    if(XpBar1(ii)<xxD1k(ii) && strcmp(geo.RaccBarrier,'ON'))
        % Calcolo approssimativo della lunghezza media di barriera:
        d1B1 = calc_distanza_punti([B1k(ii),zeros(1,nlay)],[XpBar1(ii),YpBar1(ii)]);
        d2B1 = calc_distanza_punti([xxD1k(ii),yyD1k(ii)],[XpBar1(ii),YpBar1(ii)]);
        d1B2 = calc_distanza_punti([B2k(ii),zeros(1,nlay)],[XpBar2(ii),YpBar2(ii)]);
        d2B2 = calc_distanza_punti([xxD2k(ii),yyD2k(ii)],[XpBar2(ii),YpBar2(ii)]);
        dB1=sum([d1B1,d2B1]); dB2=sum([d1B2,d2B2]); d=mean([dB1,dB2]);
        % Il raggio è variabile per le varie barriere e corrisponde ad 1/4
        % della lunghezza media di barriera...
        R_RaccB1(ii)=d/4;
        [x,y,xtan1,ytan1,xtan2,ytan2]=cir_tg_2rette(B1k(ii),0,XpBar1(ii),YpBar1(ii),XpBar1(ii),YpBar1(ii),xxD1k(ii),yyD1k(ii),R_RaccB1(ii));
        XcRacc_B1(ii)=x(1); YcRacc_B1(ii)=y(1);
        xRaccR1_B1(ii)=xtan1(1);
        yRaccR1_B1(ii)=ytan1(1);
        xRaccR2_B1(ii)=xtan2(1);
        yRaccR2_B1(ii)=ytan2(1);
        
        % Correct drawing control
        % the following script control the feasibility of the arc, if not the radius is reduced to a feasible value
        if (xRaccR2_B1(ii)>=xxD1k(ii))
            d = calc_distanza_punti([xxD1k(ii),yyD1k(ii)],[XpBar1(ii),YpBar1(ii)]);
            R_RaccB1(ii)=d/2;
            [x,y,xtan1,ytan1,xtan2,ytan2]=cir_tg_2rette(B1k(ii),0,XpBar1(ii),YpBar1(ii),XpBar1(ii),YpBar1(ii),xxD1k(ii),yyD1k(ii),R_RaccB1(ii));
            XcRacc_B1(ii)=x(1); YcRacc_B1(ii)=y(1);
            xRaccR1_B1(ii)=xtan1(1);
            yRaccR1_B1(ii)=ytan1(1);
            xRaccR2_B1(ii)=xtan2(1);
            yRaccR2_B1(ii)=ytan2(1);
            
        end
        %%
        angT2=atan2(yRaccR2_B1(ii)-YcRacc_B1(ii),xRaccR2_B1(ii)-XcRacc_B1(ii));
        angT1=atan2(yRaccR1_B1(ii)-YcRacc_B1(ii),XcRacc_B1(ii)-xRaccR1_B1(ii))+3.14;
        if (angT2<0)
            angT2=pi+angT2;
        end
        
        clear XCer YCer
        XCer=XcRacc_B1(ii)+R_RaccB1(ii)*cos(linspace(angT1,angT2,20));
        YCer=YcRacc_B1(ii)+R_RaccB1(ii)*sin(linspace(angT1,angT2,20));
    else
        R_RaccB1(ii)=0;
    end
end
%
XcRacc_B2=NaN*ones(1,nlay);
YcRacc_B2=NaN*ones(1,nlay);
xRaccR1_B2=NaN*ones(1,nlay);
yRaccR1_B2=NaN*ones(1,nlay);
xRaccR2_B2=NaN*ones(1,nlay);
yRaccR2_B2=NaN*ones(1,nlay);

for ii=1:length(B2k)
    if(XpBar2(ii)< xxD2k(ii) && strcmp(geo.RaccBarrier,'ON'))
        R_RaccB2(ii)=R_RaccB1(ii);
        [x,y,xtan1,ytan1,xtan2,ytan2]=cir_tg_2rette(B2k(ii),0,XpBar2(ii),YpBar2(ii),XpBar2(ii),YpBar2(ii),xxD2k(ii),yyD2k(ii),R_RaccB2(ii));
        % This control is necessary becouse in some cases there is an
        % intersection of the boundary of flux barrier due to the junctions,
        % statisitcally this seems to happen for the first flux barrier...
        if (ii>1 && abs((x(1)-XpBar1(ii-1))+1j*(y(1)-YpBar1(ii-1)))>R_RaccB2(ii))
            R_RaccB2(ii)=abs((x(1)-XpBar1(ii-1))+1j*(y(1)-YpBar1(ii-1)))+hfe_min;
            [x,y,xtan1,ytan1,xtan2,ytan2]=cir_tg_2rette(B2k(ii),0,XpBar2(ii),YpBar2(ii),XpBar2(ii),YpBar2(ii),xxD2k(ii),yyD2k(ii),R_RaccB2(ii));
        end
        
        XcRacc_B2(ii)=x(1); YcRacc_B2(ii)=y(1);
        xRaccR1_B2(ii)=xtan1(1);
        yRaccR1_B2(ii)=ytan1(1);
        xRaccR2_B2(ii)=xtan2(1);
        yRaccR2_B2(ii)=ytan2(1);
        
        %% Correct drawing control
        %% the following script control the feasibility of the arc, if not the radius is reduced to a feasible value
        if (xRaccR2_B2(ii)>=xxD2k(ii))
            d = calc_distanza_punti([xxD2k(ii),yyD2k(ii)],[XpBar2(ii),YpBar2(ii)]);
            R_RaccB2(ii)=d/4;
            [x,y,xtan1,ytan1,xtan2,ytan2]=cir_tg_2rette(B2k(ii),0,XpBar2(ii),YpBar2(ii),XpBar2(ii),YpBar2(ii),xxD2k(ii),yyD2k(ii),R_RaccB2(ii));
            XcRacc_B2(ii)=x(1); YcRacc_B2(ii)=y(1);
            xRaccR1_B2(ii)=xtan1(1);
            yRaccR1_B2(ii)=ytan1(1);
            xRaccR2_B2(ii)=xtan2(1);
            yRaccR2_B2(ii)=ytan2(1);
        end
        %%
        
        angT2=atan2(yRaccR2_B2(ii)-YcRacc_B2(ii),xRaccR2_B2(ii)-XcRacc_B2(ii));
        angT1=atan2(yRaccR1_B2(ii)-YcRacc_B2(ii),XcRacc_B2(ii)-xRaccR1_B2(ii))+3.14;
        if (angT2<0)
            angT2=pi+angT2;
        end
        
        clear XCer YCer
        XCer=XcRacc_B2(ii)+R_RaccB2(ii)*cos(linspace(angT1,angT2,20));
        YCer=YcRacc_B2(ii)+R_RaccB2(ii)*sin(linspace(angT1,angT2,20));
    else
        R_RaccB2(ii)=0;
    end
    
end

%%%Rotor Fillet%%%
temp.xC1k=nan(1,geo.nlay);
temp.yC1k=nan(1,geo.nlay);
temp.xC2k=nan(1,geo.nlay);
temp.yC2k=nan(1,geo.nlay);
temp.xC3k=nan(1,geo.nlay);
temp.yC3k=nan(1,geo.nlay);
temp.xC4k=nan(1,geo.nlay);
temp.yC4k=nan(1,geo.nlay);
temp.xC01k=nan(1,geo.nlay);
temp.yC01k=nan(1,geo.nlay);
temp.xC02k=nan(1,geo.nlay);
temp.yC02k=nan(1,geo.nlay);

% R466 - new rotor ends, straight + fillet 
if isfinite(geo.RotorFillet)
%>>>>>>> .r469
    if sum(RotorFillet>hc/2)
        disp('Rotor Fillet limited');
        RotorFillet(RotorFillet>hc/2)=floor(hc(RotorFillet>hc/2)/2*10^2)/10^2;
    end
    
    m1=(YpBar1-yyD1k)./(XpBar1-xxD1k);
    if ~isfinite(m1(1,1))
        m1(1,1)=m1(1,2);
    end
    q1=YpBar1-m1.*XpBar1;
    m2=-xpont./ypont;
    m1perp=-ones(1,nlay)./m1;
    m2perp=-ones(1,nlay)./m2;
    q2=ypont-m2.*xpont;
    q3=yyD2k-m1.*xxD2k;
    
    %tangential rib superior arc starts in C1k and ends in C2k with center in C01k

    xC2k=(sqrt(1+m1.*m1).*RotorFillet-q1+q2)./(m1-m2);
    yC2k=m2.*xC2k+q2;
    xC1k=(sqrt(1+m2.*m2).*RotorFillet-q2+q1)./(m2-m1);
    yC1k=m1.*xC1k+q1;
    q11=yC1k-m1perp.*xC1k;
    q22=yC2k-m2perp.*xC2k;
    xC01k=(q22-q11)./(m1perp-m2perp);
    yC01k=m2perp.*xC01k+q22;
    
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
    temp.xC1k=xC1k;
    temp.yC1k=yC1k;
    temp.xC2k=xC2k;
    temp.yC2k=yC2k;
    temp.xC3k=xC3k;
    temp.yC3k=yC3k;
    temp.xC4k=xC4k;
    temp.yC4k=yC4k;
    temp.xC01k=xC01k;
    temp.yC01k=yC01k;
    temp.xC02k=xC02k;
    temp.yC02k=yC02k;
end

% Radial ribs
temp.B1k=B1k;
temp.B2k=B2k;
temp.xpont=xpont;
temp.ypont=ypont;
temp.XpBar1=XpBar1;
temp.YpBar1=YpBar1;
temp.XpBar2=XpBar2;
temp.YpBar2=YpBar2;
temp.xxD1k=xxD1k;
temp.yyD1k=yyD1k;
temp.xxD2k=xxD2k;
temp.yyD2k=yyD2k;

geo.hc=hc;

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

temp.Bx0=Bx0;

temp.R_RaccB1=R_RaccB1; % barrier fillet (TO BE ELIMINATED)
temp.R_RaccB2=R_RaccB2; % barrier fillet
temp.XcRacc_B1=XcRacc_B1;
temp.YcRacc_B1=YcRacc_B1;
temp.xRaccR1_B1=xRaccR1_B1;
temp.yRaccR1_B1=yRaccR1_B1;
temp.xRaccR2_B1=xRaccR2_B1;
temp.yRaccR2_B1=yRaccR2_B1;
%
temp.XcRacc_B2=XcRacc_B2;
temp.YcRacc_B2=YcRacc_B2;
temp.xRaccR1_B2=xRaccR1_B2;
temp.yRaccR1_B2=yRaccR1_B2;
temp.xRaccR2_B2=xRaccR2_B2;
temp.yRaccR2_B2=yRaccR2_B2;

temp.xTraf1=xTraf1;
temp.xTraf2=xTraf2;
temp.yTraf1=yTraf1;
temp.yTraf2=yTraf2;

% temp.arcLayTraf1=arcLayTraf1;
% temp.arcLayTraf2=arcLayTraf2;

temp.XcRibTraf1=XcRibTraf1;
temp.YcRibTraf1=YcRibTraf1;
temp.XcRibTraf2=XcRibTraf2;
temp.YcRibTraf2=YcRibTraf2;

% % determination of different magnet segment, central point and magnetization direction
[geo,mat,temp] = PMdefinition_Seg(geo,mat,temp);

if strcmp(mat.LayerMag.MatName,'Air')
    temp.Mag=[];
end

% add control if the constraction of barrier is too small, so barrier is
% constructed like diamond
temp.LowDimBarrier=LowDimBarrier;
hf=[r,B1k]-[B2k,Ar]; %calcolo dei Delta fi ferro di rotore
geo.hf = hf;

% barrier transverse dimension (for permeance evaluation)
temp1_sk = calc_distanza_punti([mean([xxD1k' xxD2k'],2) mean([yyD1k' yyD2k'],2)],[mean([XpBar1' XpBar2'],2) mean([YpBar1' YpBar2'],2)]);
temp2_sk = calc_distanza_punti([mean([B1k' B2k'],2) mean([B1k' B2k'],2)*0],[mean([XpBar1' XpBar2'],2) mean([YpBar1' YpBar2'],2)]);
sk = temp1_sk' + temp2_sk' + 1/0.5 * hc/2;

geo.sk = sk;
geo.pbk = geo.sk ./ geo.hc;
geo.la = sum(geo.hc)/geo.r;
geo.lfe = sum(geo.hf)/geo.r;
geo.ly = (geo.R - (geo.r + geo.g + geo.lt))/geo.r;
geo.B1k=B1k;
geo.B2k=B2k;
geo.YpBar1=YpBar1;
geo.hc=hc;
geo.RotorFillet=RotorFillet;
