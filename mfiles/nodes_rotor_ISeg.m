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

function [geo,mat,temp]=nodes_rotor_ISeg(geo,mat)

r = geo.r;                    % Raggio del rotore al traferro
x0 = geo.x0;                    % Centro fittizio
Ar=geo.Ar;
rshaft=geo.Ar;
l = geo.l;                      % Lunghezza pacco
g = geo.g;                      % Traferro
pont0 = geo.pont0;              % minimum mechanical tolerance
pontT = geo.pontT;              % Airgap ribs [mm]

p = geo.p;                      % Paia poli
nlay = geo.nlay;                % N° layers

dalpha = geo.dalpha;            % Angoli dalpha
% Eval alpha
alpha = cumsum(dalpha);
% racc_pont = geo.racc_pont;      % racc_pont=1*pont0 <- per i ponticelli radiali.
% ang_pont0 = geo.ang_pont0;      % Ampiezza dell'angolo (in gradi) da spazzare con  raggio r in modo da ottenre un arco lungo pont0
% 
% nmax = geo.nmax;                % Velocità max (rpm) per la valutazione della sollecitazione centrifuga più gravosa (-> ponticelli)
error_mex=zeros(1,nlay);
% sigma_max = mat.Rotor.sigma_max;    % snervamento materiale [MPa]
% rhoFE = mat.Rotor.kgm3;             % densità del ferro di rotore [kg/m3]
% rhoPM = mat.LayerMag.kgm3;          % densità magneti [kg/m3]

%% Determination of air thickness and check the feasibility of the geometry
geo = calcHcCheckGeoControl(geo);
hc=geo.hc;

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

B1k_temp=geo.B1k;
B2k_temp=geo.B2k;
hfe=[r-B2k_temp(1),B1k_temp(1:nlay-1)-B2k_temp(2:nlay),B1k_temp(nlay)-Ar]; %%

r_all = zeros(2*nlay,1);                                                                 % Nel vettore r_all si memorizzano invece le distanze, sempre rispetto al
% Starting and ending pooint of the flux barrier...
for jj = 1:nlay                                                             % che individuano l'inizio e la fine delle nlay barriere.
    r_all([2*jj-1 2*jj]) = [x0-B2k_temp(jj) x0-B1k_temp(jj)];
end
% 
% pont_bound_conflict1=find(r_all(1:2:2*nlay)>=rbeta);
% pont_bound_conflict2=find(r_all(2:2:2*nlay)<=rbeta);
%% Posizione banane su asse q (convenzione assi VAGATI)
% 2013/07/06 MG punto banana su asse q serve per i successivi export in DXF

% hcc=(r_all(2:2:end))-(r_all(1:2:end));
[xpont,ypont] = calc_intersezione_cerchi((r-pontT), rbeta, x0);    % valore non del tutto corretto ricalcolato in seguito, tuttavia serve per il disegno della curva centro bar...
[xcbar,ycbar] = calc_intersezione_cerchi((r-pontT-hc/2), rbeta, x0);
% 2014/02/26 MG Problem of immaginary number for higher nlay value
for ii=1:nlay
    if (not(isreal(xpont(ii)))||not(isreal(ypont(ii))))
        xpont(ii)=r-pontT(ii); ypont(ii)=0;
        error_mex(ii)=1;
    end
    if (not(isreal(xcbar(ii)))||not(isreal(ycbar(ii))))
        xcbar(ii)=r-pontT(ii)-hc(ii)/2;
        ycbar(ii)=0;
        error_mex(ii)=1;
    end
end

rpont_x0=sqrt(ypont.^2+(x0-xpont).^2);
Bx0=x0-(rpont_x0);

%% Calc di xpont, ypont:
for ii=1:nlay
    [a,b,c]=retta_per_2pti(0,0,xcbar(ii),ycbar(ii));
    if (a==0 && c==0)
        % xpont and ypont are just evaluated
    else
        A=1+b^2/a^2; B=2*b*c/a; C=(c^2/a^2-(r-pontT(ii))^2);
        ytemp=roots([A,B,C]); ypont(ii)=ytemp(find(ytemp>=0));
        xpont(ii)=-(b*ypont(ii)+c)/a;
    end
end

for ii=1:nlay
    if ii==1
        B2k(ii)=xcbar(ii)+hc(ii)/2;
        B1k(ii)=xcbar(ii)-hc(ii)/2;        
    else
        B1k(ii)=B1k(ii-1)-hfe(ii)-hc(ii);
        B2k(ii)=B1k(ii-1)-hfe(ii);
    end
    
end

%% Intersezione circonferenze punti al traferro:
[xTraf2,yTraf2] = calc_intersezione_cerchi(r-pontT, r_all(1:2:end), x0);
[xTraf1,yTraf1] = calc_intersezione_cerchi(r-pontT, r_all(2:2:end), x0);

for ii=1:nlay
    if (not(isreal(yTraf2(ii))))
        xcbar(ii)=r-pontT(ii)-hc(ii)/2; ycbar(ii)=0;
        error_mex(ii)=1;
    end
end
%% 2014/06/22 MG
% set di istruzioni per il disegno delle barriere nel caso in cui la punta
% della barriera sia interna alla ipotetica curva centrale di barriera...
for ii=2:nlay
    if (xTraf2(ii)<=xpont(ii) || yTraf2(ii)>=ypont(ii))
        [xpont_temp,ypont_temp] = calc_intersezione_cerchi((r-pontT(ii)), rbeta(ii), x0);
        xpont(ii)=xpont_temp;
        ypont(ii)=ypont_temp;
        clear xpont_temp ypont_temp;
    end
end
%% Intersezione tra rette che compongono i lati delle barriere di flux:
% retta verticale
for ii=1:nlay
    if ii==1
        XpBar2(ii)=xpont(1);
        YpBar2(ii)=ypont(1);
        
    else
        a1=1;
        b1=0;
        c1=-B2k(ii);
        m2=tan(pi/2/p);
        a2=m2;
        b2=-1;
        c2=(yTraf2(ii)-m2*xTraf2(ii));
        [x,y]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
        XpBar2(ii)=x;
        YpBar2(ii)=y;
    end
end

for ii=1:nlay
    if ii==1
        XpBar1(ii)=B1k(1);
        YpBar1(ii)=ypont(1);
        
    else
        a1=1;
        b1=0;
        c1=-B1k(ii);
        m2=tan(pi/2/p);
        a2=m2;
        b2=-1;
        c2=(yTraf1(ii)-m2*xTraf1(ii));
        [x,y]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
        XpBar1(ii)=x;
        YpBar1(ii)=y;
    end
end


for ii=1:nlay
    
    if ii==1
        
        xD1k=B1k(1);
        yD1k=ycbar(1);
        xxD1k(ii)=B1k(1);
        yyD1k(ii)=ycbar(1);
        XpBar1(1)=xD1k;
        YpBar1(1)=yD1k;
        
        rcir=hc(1)/2;
        xcir=xcbar(1); ycir=ycbar(1);
        
        angT2=atan2(yD1k-ycir,xD1k-xcir);
        angT1=atan2(hc(1)/2,0);
        if (angT2<0)
            angT2=pi+angT2;
        end
        
        arcLayTraf1(ii)=(angT2-angT1)*180/pi;
        XcRibTraf1(ii)=xcbar(1);
        YcRibTraf1(ii)=ycbar(1);
        
        clear XCer YCer
        %             XCer=xcir+rcir*cos(linspace(angT1,angT2,20));
        %             YCer=ycir+rcir*sin(linspace(angT1,angT2,20));
        %             figure(100);hold on;plot(XCer,YCer,'+m');hold off;
        %             figure(100);hold on;plot(xD1k,yD1k,'rs');plot(xcbar(ii),ycbar(ii),'gs'); hold off;
        %             figure(100);hold on;plot([0,xpont(ii)]',[0,ypont(ii)]','-b'); hold off;
        
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
        %             XCer=xc+rc*cos(linspace(angT1,angT2,20));
        %             YCer=yc+rc*sin(linspace(angT1,angT2,20));
        %             figure(100);hold on;plot(XCer,YCer,'+m');hold off;
        %             figure(100);hold on;plot(xD1k,yD1k,'rs');plot(xc,yc,'gs'); hold off;
        %             figure(100);hold on;plot([0,xpont(ii)]',[0,ypont(ii)]','-b'); hold off;
    end
end

for ii=1:nlay
    if ii==1
        xD2k=B2k(1);
        yD2k=ycbar(1);
        xxD2k(ii)=xD2k;
        yyD2k(ii)=yD2k;
        XpBar2(1)=xD2k;
        YpBar2(1)=yD2k;
        
        rcir=hc(1)/2;
        xcir=xcbar(1); ycir=ycbar(1);
        
        angT2=atan2(yD2k-ycir,xD2k-xcir);
        angT1=atan2(hc(1)/2,0);
        if (angT2<0)
            angT2=pi+angT2;
        end
        
        arcLayTraf2(ii)=(angT2-angT1)*180/pi;
        XcRibTraf2(ii)=xcbar(1);
        YcRibTraf2(ii)=ycbar(1);
        
        clear XCer YCer
        %             XCer=xcir+rcir*cos(linspace(angT1,angT2,20));
        %             YCer=ycir+rcir*sin(linspace(angT1,angT2,20));
        %             figure(100);hold on;plot(XCer,YCer,'+m');hold off;
        %             figure(100);hold on;plot(xD1k,yD1k,'rs');plot(xcbar(ii),ycbar(ii),'gs'); hold off;
        %             figure(100);hold on;plot([0,xpont(ii)]',[0,ypont(ii)]','-b'); hold off;
        %
        
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
        %             XCer=xc+rc*cos(linspace(angT1,angT2,20));
        %             YCer=yc+rc*sin(linspace(angT1,angT2,20));
        %             figure(100);hold on;plot(XCer,YCer,'+m');hold off;
    end
    
end
%%
%% Controllo di sicurezza per evitare che la circonferenza del rib traf termini dopo il lato di barriera
%%
for ii=2:nlay
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
    
    if (yyD1k(ii)<=YpBar1(ii) && xxD1k(ii)<=XpBar1(ii))
        yyD1k(ii)=YcRibTraf1(ii)+sqrt((abs((xpont(ii)-XcRibTraf1(ii))+1j*(ypont(ii)-YcRibTraf1(ii))))^2-(XpBar1(ii)-XcRibTraf1(ii))^2);
        YpBar1(ii)=yyD1k(ii);
        xxD1k(ii)=XpBar1(ii);
    end
    
end

if (error_mex(1)==1)
    temp.xc=(xxD1k(1:end)+xxD2k(1:end))/2;
    temp.yc=[0,(yyD1k(2:end)+yyD2k(2:end))/2];
else
    temp.xc=(xxD1k+xxD2k)/2;
    temp.yc=(yyD1k+yyD2k)/2;
end

%%
%% Disegno del lamierino cn figure...
%%
% %%
% xcir_plot=[0:0.5:r];
% ycir_plot=sqrt(r^2-xcir_plot.^2);
% figure(100);hold on;
% xo=x0-rpont_x0'*cos(0:0.1:pi/2);
% yo=rpont_x0'*sin(0:0.1:pi/2);
% plot(xo',yo','--r','LineWidth',2); axis([0 r+0.5 0 r+0.5]); axis square
% plot(B1k,0,'ob');plot(B2k,0,'ob');
% plot(xpont,ypont,'*c');
% plot(XpBar2,YpBar2,'bs');
% plot(XpBar1,YpBar1,'bs');
% plot(xTraf1,yTraf1,'*m');
% plot(xTraf2,yTraf2,'*m');
% plot(xD1k,yD1k,'r^');
% plot(xcir_plot,ycir_plot,'k');
%
% for ii=1:nlay
%    plot([XpBar1(ii),xxD1k(ii)],[YpBar1(ii),yyD1k(ii)],'b','LineWidth',2);
%    plot([XpBar2(ii),xxD2k(ii)],[YpBar2(ii),yyD2k(ii)],'b','LineWidth',2);
%    plot([B1k(ii),XpBar1(ii)],[0,YpBar1(ii)],'b','LineWidth',2);
%    plot([B2k(ii),XpBar2(ii)],[0,YpBar2(ii)],'b','LineWidth',2);
%
% end
% hold off;
% keyboard
%
%
%
% Determinazione dei punti caratteristici per il calcolo delle perdite nel
% ferro:
calcolo_posizioni_Pfe;
% radial ribs evaluation:
rTemp = rbeta;

temp.Bx0=Bx0;

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

[temp,geo] = calc_ribs_rad_Seg(geo,mat,temp);

%calc_ribs_rad;

% Additional division for magnet insertion in flux barrier
% XpMag1B1=XpontRadBarSx;
% YpMag1B1=YpBar2;
% temp.XpMag1B1=XpMag1B1;
% temp.YpMag1B1=YpMag1B1;

% output
% temp.B1k=B1k;
% temp.B2k=B2k;
% temp.Bx0=Bx0;
% temp.xpont=xpont;
% temp.ypont=ypont;
% temp.xxD1k=xxD1k;
% temp.yyD1k=yyD1k;
% temp.xxD2k=xxD2k;
% temp.yyD2k=yyD2k;
% temp.XpBar1=XpBar1;
% temp.YpBar1=YpBar1;
% temp.XpBar2=XpBar2;
% temp.YpBar2=YpBar2;

temp.xTraf1=xTraf1;
temp.xTraf2=xTraf2;
temp.yTraf1=yTraf1;
temp.yTraf2=yTraf2;

temp.arcLayTraf1=arcLayTraf1;
temp.arcLayTraf2=arcLayTraf2;

temp.XcRibTraf1=XcRibTraf1;
temp.YcRibTraf1=YcRibTraf1;
temp.XcRibTraf2=XcRibTraf2;
temp.YcRibTraf2=YcRibTraf2;

% Points for radial ribs
% temp.XpontRadDx=XpontRadDx;
% temp.YpontRadDx=YpontRadDx;
% temp.XpontRadSx=XpontRadSx;
% temp.YpontRadSx=YpontRadSx;
% temp.XpontRadBarDx=XpontRadBarDx;
% temp.XpontRadBarSx=XpontRadBarSx;
% temp.YpontRadBarDx=YpontRadBarDx;
% temp.YpontRadBarSx=YpontRadBarSx;
temp.error_mex=error_mex;
%Aree e matrice Mag per posizionamento magneti

[geo,mat,temp] = PMdefinition_Seg(geo,mat,temp);

if strcmp(mat.LayerMag.MatName,'Air')
    temp.Mag=[];
end

% [temp,geo]=area_magnet_ISeg(temp,geo);
% % determination of different magnet segment, central point and magnetization direction
% MagnetFullFill_ISeg;
% 
% mat.LayerMag.Br = [Br Br];   % doubles Br pieces (half pole + half pole)
% 
% temp.xc=xc;
% temp.yc=yc;
% temp.xair=xair;
% temp.yair=yair;
% temp.xmag=xmag;
% temp.ymag=ymag;
% temp.xmagair=xmagair;
% temp.ymagair=ymagair;
% temp.zmag=zmag;
% beta_f=0;
hf=[r,B1k]-[B2k,Ar]; % calcolo dei Delta di ferro di rotore
geo.hf = hf;
% geo.beta_f = beta_f;
%geo.pontR=pont;

% barrier transverse dimension (for permeance evaluation)
temp1_sk = calc_distanza_punti([mean([xxD1k' xxD2k'],2) mean([yyD1k' yyD2k'],2)],[mean([XpBar1' XpBar2'],2) mean([YpBar1' YpBar2'],2)]);
temp2_sk = calc_distanza_punti([mean([B1k' B2k'],2) mean([B1k' B2k'],2)*0],[mean([XpBar1' XpBar2'],2) mean([YpBar1' YpBar2'],2)]);
% equivalent sk = segment 1 + segment 2 + barrier end radius weighted by 1/0.7822
sk = temp1_sk'+temp2_sk' + 1/0.7822 * hc/2;

geo.sk = sk;
geo.pbk = geo.sk ./ geo.hc;
geo.la = sum(geo.hc)/geo.r;
geo.lfe = sum(geo.hf)/geo.r;
geo.ly = (geo.R - (geo.r + geo.g + geo.lt))/geo.r;


end

