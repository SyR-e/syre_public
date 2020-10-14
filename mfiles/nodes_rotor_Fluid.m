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

function [geo,mat,temp]=nodes_rotor_Fluid(geo,mat)

l=geo.l;
r = geo.r;            % Raggio del rotore al traferro
rlim=r;
p = geo.p;              % Paia poli
nlay = geo.nlay;        % N° layers
R = geo.R;              % Raggio ext
g = geo.g;              % Traferro
lt = geo.lt;            % Lunghezza denti
pont0 = geo.pont0;      % Ponticelli al traferro
rshaft=geo.Ar;          % raggio albero fittizio per le linee di lfux
Ar=geo.Ar;              % raggio albero
nmax=geo.nmax;          % velocità massima
BandRib=[2.6*geo.pont0,5*geo.pont0,5*geo.pont0];
rCirRib=1*pont0;
sigma_max = mat.Rotor.sigma_max;    % snervamento materiale [MPa]
rhoFE = mat.Rotor.kgm3;             % densità del ferro di rotore [kg/m3]
rhoPM = mat.LayerMag.kgm3;          % densità magneti [kg/m3]

x0 = r/cos(pi/2/p);                            % centro cerchi (p generico)
geo.x0 = x0;
dx=geo.dx;
dalpha = geo.dalpha;            % Angoli dalpha
% Eval alpha
alpha = cumsum(dalpha);
alphaRot=alpha*pi/180+pi/2/p;   % angoli barriere di flux al traferro meccanici
LastBarCurvatura=1;
%%
%% Determinazione delle costanti per il disegno delle linee mediane in
%%
% corrispondenza degli alphaRot:
C=sin(p*alphaRot).*((((r-pont0)./rshaft).^(2*p)-1)./((r-pont0)./rshaft).^p);

Dteta=1;
teta=[0:Dteta*pi/180:pi/p];
rTemp=zeros(length(C),length(teta));
z=zeros(length(C),length(teta));
B0=zeros(length(C),1);

for k=1:length(C)
    
    rTemp(k,:)=rshaft*((C(k)+sqrt(C(k)^2+4*(sin(p*teta)).^2))./(2*sin(p*teta))).^(1/p);
    z(k,:)=rTemp(k,:).*exp(1j*teta);
    Br0(k,:)=rshaft*((C(k)+sqrt(C(k)^2+4))./2).^(1/p);
    % coordinata angolare a partire dall'asse y centrale meno il ponticello.
    tetaRpont(k)=180/p-(1/p)*asin((C(k)*((r-geo.pont0)/rshaft)^p)/(((r-geo.pont0)/rshaft)^(2*p)-1))*180/pi;
    tetaRpont1(k)=180/p-(1/p)*asin((C(k)*((r-geo.pont0-rCirRib)/rshaft)^p)/(((r-geo.pont0-rCirRib)/rshaft)^(2*p)-1))*180/pi;
    
end

x=real(z); y=imag(z);
[xo,yo]=rot_point(x',y',-pi/2/p);

% coordinate in x,y al traferro meno il ponticello
[xpont,ypont]=rot_point((r-geo.pont0)*cos(tetaRpont*pi/180),(r-geo.pont0)*sin(tetaRpont*pi/180),-pi/2/p);

%% rotazione di -pi/2/p (rotazione nel 1° e 4° quadrante):
[Bx0,By0]=rot_point(Br0*cos(pi/2/p),Br0*sin(pi/2/p),-pi/2/p);
Bx0=Bx0';
By0=By0';
%%
%% Determination of air thickness and check the feasibility of the geometry
geo.Bx0=Bx0; % Initialization of central non-moved line of the flux barrier
geo = calcHcCheckGeoControlwDx(geo);
%%
%  hc=geo.hc;
B1k=geo.B1k;
B2k=geo.B2k;
hc=B1k-B2k;
% Instruction for flux barrier translaction in function of range dx=[-1:1]:
for k=1:length(dx)
    if (dx(k)==1)||(dx(k)==-1)
        Bx0(k)=(B1k(k)+B2k(k))/2;
        [xBx0New(k),yBy0New(k)]=rot_point(Bx0(k),0,pi/2/p);
        Br0New(k)=abs(xBx0New(k)+1j*yBy0New(k));
        CNew(k)=(((Br0New(k)./rshaft).^(2*p)-1)./(Br0New(k)./rshaft).^p);
        tetaRpontNew(k)=180/p-(1/p)*asin((CNew(k)*((r-geo.pont0)/rshaft)^p)/(((r-geo.pont0)/rshaft)^(2*p)-1))*180/pi;
        % Riassegno xpont ed ypont
        [xpontNew,ypontNew]=rot_point((r-geo.pont0)*cos(tetaRpontNew(k)*pi/180),(r-geo.pont0)*sin(tetaRpontNew(k)*pi/180),-pi/2/p);
        xpont(k)=xpontNew;
        ypont(k)=ypontNew;
    end
end

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%
%% SECURITY CONTROL END
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%

%% ROTOR LAYER DRAWING:
% si ruotano i valori di Bk nel 1° quadrante ( poichè la validità delle formule per la forma delle barriere è limitata al 1° quadrante):

[xB1k,yB1k]=rot_point(B1k,0,pi/2/p);
[xB2k,yB2k]=rot_point(B2k,0,pi/2/p);

rB1k=abs(xB1k+1j*yB1k);
rB2k=abs(xB2k+1j*yB2k);
CB1k=(((rB1k./rshaft).^(2*p)-1)./(rB1k./rshaft).^p);
CB2k=(((rB2k./rshaft).^(2*p)-1)./(rB2k./rshaft).^p);

% Ciclo per il disegno dei lati delle barriere di flux:
for k=1:length(CB1k)
    
    RB1k(k,:)=rshaft*((CB1k(k)+sqrt(CB1k(k)^2+4*(sin(p*teta)).^2))./(2*sin(p*teta))).^(1/p);
    zB1k(k,:)=RB1k(k,:).*exp(1j*teta);
    
    RB2k(k,:)=rshaft*((CB2k(k)+sqrt(CB2k(k)^2+4*(sin(p*teta)).^2))./(2*sin(p*teta))).^(1/p);
    zB2k(k,:)=RB2k(k,:).*exp(1j*teta);
    
    tetaIBk1(k)=180/p-(1/p)*asin((CB1k(k)*((r-geo.pont0)/rshaft)^p)/(((r-geo.pont0)/rshaft)^(2*p)-1))*180/pi;
    tetaIBk2(k)=180/p-(1/p)*asin((CB2k(k)*((r-geo.pont0)/rshaft)^p)/(((r-geo.pont0)/rshaft)^(2*p)-1))*180/pi;
    tetaIBK2=45;
    tetaTraf01(k)=180/p-(1/p)*asin((CB1k(k)*(r/rshaft)^p)/((r/rshaft)^(2*p)-1))*180/pi;
    tetaTraf02(k)=180/p-(1/p)*asin((CB2k(k)*(r/rshaft)^p)/((r/rshaft)^(2*p)-1))*180/pi;
    
end

% rotazione dei valori nel 1°e 4° quadrante pronti per il disegno in matlab
% del lamierino (no FEMM)

[xxB1k,yyB1k]=rot_point(real(zB1k(:,2:end)),imag(zB1k(:,2:end)),-pi/2/p);
[xxB2k,yyB2k]=rot_point(real(zB2k(:,2:end)),imag(zB2k(:,2:end)),-pi/2/p);

% Determinazione xTraf e yTraf punti al traferro:
[xTraf01,yTraf01]=rot_point(r*cos(tetaTraf01*pi/180),r*sin(tetaTraf01*pi/180),-pi/2/p);
[xTraf02,yTraf02]=rot_point(r*cos(tetaTraf02*pi/180),r*sin(tetaTraf02*pi/180),-pi/2/p);
% Determinazione xTraf e yTraf punti al traferro - pont0:
[xTraf1,yTraf1]=rot_point((r-geo.pont0)*cos(tetaIBk1*pi/180),(r-geo.pont0)*sin(tetaIBk1*pi/180),-pi/2/p);
[xTraf2,yTraf2]=rot_point((r-geo.pont0)*cos(tetaIBk2*pi/180),(r-geo.pont0)*sin(tetaIBk2*pi/180),-pi/2/p);

%% CHECK of low flux barrier length
LowBarLength=0;
if (yTraf1(1)<=eps)
    yTraf1(1)=eps;
    LowBarLength=1;
elseif (yTraf2(1)<=eps)
    yTraf2(1)=eps;
    LowBarLength=1;
end
%% RE-Assignation and interpolation of flux boundary
% per la parte positiva in xy di geo mot
[x,y]=interp_flux_barrier(xxB1k,yyB1k,2*yTraf01);
clear xxB1k yyB1k;
xxB1k=x; yyB1k=y;
clear x y
[x,y]=interp_flux_barrier(xxB2k,yyB2k,2*yTraf02);
clear xxB2k yyB2k;
xxB2k=x; yyB2k=y;

%%
%% DISEGNO DEI PUNTI AL TRAFERRO MEDIANTE TANGENTE TRA CIRCONFERENZE:
%%

for k=1:length(CB2k)
    if (k==1 && LowBarLength==1)
        xxD2k(1)=xpont(1);
        yyD2k(1)=ypont(1);
    else
        [xxB2k_med,yyB2k_med]=valore_medio_di_barriera(xxB2k(k,:),yyB2k(k,:),yTraf2(k));
        [x02,y02,r02]=circonferenza_per_3_pti(B2k(k),0,xTraf2(k),yTraf2(k),xxB2k_med,yyB2k_med);
        rTemp=r-pont0;
        [xD2k,yD2k,xc,yc,rc]=cir_tg_2cir(xpont(k),ypont(k),rTemp,x02,y02,r02);
        angT1=atan2(yD2k-yc,xD2k-xc);
        angT2=atan2(ypont(k)-yc,xpont(k)-xc);
        %     if (angT1<angT2)
        %         angT2=-2*pi+angT2;
        %     end
        XcRibTraf2(k)=xc;
        YcRibTraf2(k)=yc;
        arcLayTraf2(k)=(angT2-angT1)*180/pi;
        xxD2k(k)=xD2k;
        yyD2k(k)=yD2k;
        clear XCer YCer
        XCer=xc+rc*cos([angT1:0.01:angT2]);
        YCer=yc+rc*sin([angT1:0.01:angT2]);
        %     figure(100);plot(XCer,YCer,'k');hold on;
    end
end
%% 2014/08/06 MG check if xxD2k(1)<xpont to avoid too low circle
if xxD2k(1)<xpont(1)
    xxD2k(1)=xpont(1);
end

for k=1:length(CB1k)
    if (k==1 && LowBarLength==1)
        xxD2k(1)=xpont(1);
        yyD2k(1)=ypont(1);
    else
        [xxB1k_med,yyB1k_med]=valore_medio_di_barriera(xxB1k(k,:),yyB1k(k,:),yTraf1(k));
        [x02,y02,r02,angleA,angleB]=circonferenza_per_3_pti(B1k(k),0,xTraf1(k),yTraf1(k),xxB1k_med,yyB1k_med);
        if (angleA<angleB)
            angleB=-2*pi+angleB;
        end
        XCerchio=x02+r02*cos([angleB:0.01:angleA]);
        YCerchio=y02+r02*sin([angleB:0.01:angleA]);
        %     figure(100);hold on;plot(XCerchio,YCerchio,'k');hold off;
        
        rTemp=r-pont0;
        [xD1k,yD1k,xc,yc,rc]=cir_tg_2cir(xpont(k),ypont(k),rTemp,x02,y02,r02);
        angT2=atan2(yD1k-yc,xD1k-xc);
        angT1=atan2(ypont(k)-yc,xpont(k)-xc);
        %     figure(100);hold on;plot(xD1k,yD1k,'rs');plot(xc,yc,'gs'); hold off;
        %     figure(100);hold on;plot([0,xpont(k)]',[0,ypont(k)]','-b'); hold off;
        
        if (angT2<0)
            angT2=pi+angT2;
        end
        arcLayTraf1(k)=(angT2-angT1)*180/pi;
        XcRibTraf1(k)=xc;
        YcRibTraf1(k)=yc;
        xxD1k(k)=xD1k;
        yyD1k(k)=yD1k;
        clear XCer YCer
        XCer=xc+rc*cos(linspace(angT1,angT2,20));
        YCer=yc+rc*sin(linspace(angT1,angT2,20));
        %     figure(100);hold on;plot(XCer,YCer,'k');hold off;
    end
end

clear xxB1k_mean yyB1k_mean xxB2k_mean yyB2k_mean;

% figure(100);hold on;
% plot(xo,yo,'--r','LineWidth',2); axis([0 rlim 0 rlim]); axis square
% % plot(r*cos(geo.alpha*pi/180),r*sin(geo.alpha*pi/180),'or');
% plot(xpont,ypont,'or');
% plot(xTraf1,yTraf1,'ob',xTraf2,yTraf2,'ob');
% plot(xxD1k,yyD1k,'bs');
% plot(xxD2k,yyD2k,'bs');
% xx=(0:0.1:rlim);
% plot(xx,sqrt(rlim^2-xx.^2),'k','LineWidth',2);
% plot([0:0.1:rshaft],sqrt(rshaft^2-[0:0.1:rshaft].^2),'k','LineWidth',2);
% plot(xxB1k',yyB1k','b','LineWidth',2);plot(xxB2k',yyB2k','b','LineWidth',2);
% plot(B1k,0,'ob');plot(B2k,0,'ob');
% % plot(xRib1,yRib1,'ms',xRib2,yRib2,'ms');
% % plot(xpont1,ypont1,'om');
% hold off;grid on;
%%
%% SET DI ISTRUZIONI INTERPOLAZIONE BARRIERE DI FLUX:
%%
[xxB1k_mean,yyB1k_mean]=valore_medio_di_barriera(xxB1k,yyB1k,yTraf1);
[xxB2k_mean,yyB2k_mean]=valore_medio_di_barriera(xxB2k,yyB2k,yTraf2);
% figure(100);hold on; plot(xxB1k_mean,yyB1k_mean,'og',xxB2k_mean,yyB2k_mean,'om');hold off;
%%

for i=1:nlay
    if (i==1)
        if (yyB2k(i,2)>=yyD2k(i)||yyD2k(1)<=0||isnan(yyD2k(1)))
            tetaBuccia2k=180/p-(1/p)*asin((CB2k(i)*((r-geo.pont0-BandRib(i))/rshaft)^p)/(((r-geo.pont0-BandRib(i))/rshaft)^(2*p)-1))*180/pi;
            % rotazione nel 1° e 4° quadrante:
            [xBuccia2k,yBuccia2k]=rot_point((r-geo.pont0-BandRib(i))*cos(tetaBuccia2k*pi/180),(r-geo.pont0-BandRib(i))*sin(tetaBuccia2k*pi/180),-pi/2/p);
            xxD2k(1)=xBuccia2k;
            yyD2k(1)=yBuccia2k;
        else
            [xxB2k_mean2,yyB2k_mean2]=valore_medio_di_barriera(xxB2k(i,:),yyB2k(i,:),yyD2k(i));
            [xc,yc,rTemp,angleA,angleB]=circonferenza_per_3_pti(B2k(i),0,xxD2k(i),yyD2k(i),xxB2k_mean2,yyB2k_mean2);
            % 2014/05/02 MG se l'archetto rotorico è troppo piccolo il
            % sistema di interpolazione sbaglia e di seguito riporto una
            % soluzione tampone per il momento:
            if (xc<=B2k(i))
                [xc,yc,rTemp,angleA,angleB]=circonferenza_per_3_pti(B2k(i),0,xxD2k(i),yyD2k(i),xTraf2(i),yTraf2(i));
            end
            
            if (angleA<angleB)
                angleB=-2*pi+angleB;
            end
            XCerchio=xc+rTemp*cos(linspace(angleB,angleA,20));
            YCerchio=yc+rTemp*sin(linspace(angleB,angleA,20));
            XcBar2(i)=xc;
            YcBar2(i)=yc;
            %             figure(100);hold on;plot(XCerchio,YCerchio,'*k');hold off;
            arcLayer2(i)=(angleA-angleB)*180/pi;
            if abs(arcLayer2(i))>90
                arcLayer2(i)=360-abs(arcLayer2(i));
            end
            
        end
    else
        [xc,yc,rTemp,angleA,angleB]=circonferenza_per_3_pti(B2k(i),0,xxD2k(i),yyD2k(i),xxB2k_mean(i),yyB2k_mean(i));
        
        if (angleA<angleB)
            angleB=-2*pi+angleB;
        end
        XCerchio=xc+rTemp*cos(linspace(angleB,angleA,20));
        YCerchio=yc+rTemp*sin(linspace(angleB,angleA,20));
        XcBar2(i)=xc;
        YcBar2(i)=yc;
        %         figure(100);hold on;plot(XCerchio,YCerchio,'k');hold off;
        arcLayer2(i)=(angleA-angleB)*180/pi;
        if abs(arcLayer2(i))>90
            arcLayer2(i)=360-abs(arcLayer2(i));
        end
        
    end % (i==1 && LowBarLength==0)
    
end % for

if nlay~=1
    
    for i=1:nlay-1
        [xc,yc,rTemp,angleA,angleB]=circonferenza_per_3_pti(B1k(i),0,xxD1k(i),yyD1k(i),xxB1k_mean(i),yyB1k_mean(i));
        
        if (angleA<angleB)
            angleB=-2*pi+angleB;
        end
        
        XCerchio=xc+rTemp*cos(linspace(angleB,angleA,20));
        YCerchio=yc+rTemp*sin(linspace(angleB,angleA,20));
        
        XcBar1(i)=xc;
        YcBar1(i)=yc;
        
        %     figure(100);hold on;plot(XCerchio,YCerchio,'k');hold off;
        arcLayer1(i)=(angleA-angleB)*180/pi;
    end
    %%
    %% SET ISTRUZIONI INTERPOLAZIONI ULTIMA BARRIERA DI FLUX:
    %%
    [xxB1k_mean2,yyB1k_mean2]=valore_medio_di_barriera(xxB1k(nlay,:),yyB1k(nlay,:),yyB1k_mean(end));
    
    % figure(100);hold on; plot(xxB1k_mean2,yyB1k_mean2,'bs');hold off;
    clear index;
    
    %% 2013/10/04 MG determinazione del segmento intermedio
    
    index=find(yyB1k(nlay,:)>=(yyB1k_mean(end)) & yyB1k(nlay,:)<=0.9*yyD1k(end));
    
    if length(index)==1
        index=find(yyB1k(nlay,:)>=(yyB1k_mean(end)) & yyB1k(nlay,:)<=1.1*yyD1k(end));
    end
    
    yyB1kAux3=yyB1k(nlay,index);
    xxB1kAux3=xxB1k(nlay,index);
    [xxB1k_mean3,yyB1k_mean3]=valore_medio_di_barriera(xxB1kAux3,yyB1kAux3,0.9*yyD1k(end));
    
    % figure(100);hold on; plot(xxB1k_mean3,yyB1k_mean3,'k');hold off;
    clear index;
    
    if LastBarCurvatura==1
        [xc,yc,rTemp,angleA,angleB]=circonferenza_per_3_pti(xxB1k_mean2,yyB1k_mean2,xxD1k(end),yyD1k(end),xxB1k_mean3,yyB1k_mean3);
    else
        % [xc,yc,rTemp,angleA,angleB]=circonferenza_per_3_pti(xxB1k_mean2,yyB1k_mean2,xxB1k_mean(end),yyB1k_mean(end),xxB1k_mean3,yyB1k_mean3);
        [xc,yc,rTemp,angleA,angleB]=circonferenza_per_3_pti(xxB1k_mean2,yyB1k_mean2,xxB1k_mean3,yyB1k_mean3,xxB1k_mean(end),yyB1k_mean(end));
    end
    % centro raggio ultima barriera di flux:
    XcBarLast_mean=xc;
    YcBarLast_mean=yc;
    
    if (angleA<angleB)
        angleB=-2*pi+angleB;
    end
    XCerchio=xc+rTemp*cos([ angleB:0.01:angleA]);
    YCerchio=yc+rTemp*sin([ angleB:0.01:angleA]);
    
    %     figure(100);hold on;plot(XCerchio,YCerchio,'*k');hold off;
    arcLayer3=(angleA-angleB)*180/pi;
    [xD1k,yD1k,xc,yc,rc]=tg_cir(xxB1k_mean3,yyB1k_mean3,xTraf1(end),yTraf1(end),xpont(end),ypont(end));
    XcRibTraf1(end)=xc;
    YcRibTraf1(end)=yc;
    angT2=atan2(yD1k-yc,xD1k-xc);
    angT1=atan2(ypont(k)-yc,xpont(k)-xc);
    %     figure(100);hold on;plot(xD1k,yD1k,'k');plot(xc,yc,'k'); hold off;
    %     figure(100);hold on;plot([0,xpont(k)]',[0,ypont(k)]','k'); hold off;
    
    if (angT2<0)
        angT2=pi+angT2;
    end
    arcLayTraf1(end)=(angT2-angT1)*180/pi;
    xxD1k(end)=xD1k;
    yyD1k(end)=yD1k;
    clear XCer YCer
    XCer=xc+rc*cos(linspace(angT1,angT2,20));
    YCer=yc+rc*sin(linspace(angT1,angT2,20));
    %     figure(100);hold on;plot(XCer,YCer,'k');hold off;
    
else
    [xc,yc,rTemp,angleA,angleB]=circonferenza_per_3_pti(B1k(1),0,xxD1k(1),yyD1k(1),xxB1k_mean(1),yyB1k_mean(1));
    
    if (angleA<angleB)
        angleB=-2*pi+angleB;
    end
    
    XCerchio=xc+rTemp*cos(linspace(angleB,angleA,20));
    YCerchio=yc+rTemp*sin(linspace(angleB,angleA,20));
    
    XcBar1(1)=xc;
    YcBar1(1)=yc;
    
    %     figure(100);hold on;plot(XCerchio,YCerchio,'k');hold off;
    arcLayer1(1)=(angleA-angleB)*180/pi;
    
    xxB1k_mean2=NaN;
    yyB1k_mean2=NaN;
    xxB1k_mean3=NaN;
    yyB1k_mean3=NaN;
    XcBarLast_mean=NaN;
    YcBarLast_mean=NaN;
    arcLayer3=NaN;
    
end

%%
%%
%% Calcolo dei volumi di ferro:
%%
calc_Iron_area_fluid;
%%
% Determinazione dei punti caratteristici per il calcolo delle perdite nel
% ferro:
% calcolo_posizioni_Pfe;

%% Assigment of central flux barrier coordinate for label
xcbar=(xxD1k+xxD2k)/2;
ycbar=(yyD1k+yyD2k)/2;
% if (yTraf1(1)<=eps || yTraf2(1)<=eps||xxD2k(1)==xpont(1)||angle((xpont(1)-xxD1k(1))+1j*(ypont(1)-yyD1k(1)))>80*pi/180 || LowBarLength==1)
%     xcbar=(xxD1k(2:end)+xxD2k(2:end))/2;
%     ycbar=(yyD1k(2:end)+yyD2k(2:end))/2;
% else
%     xcbar=(xxD1k+xxD2k)/2;
%     ycbar=(yyD1k+yyD2k)/2;
% end
% determination of different magnet segment, central point and
% magnetization direction
% geo.Br=0.4;
xmag=[];
ymag=[]; Br = [];

for kk=1:length(xcbar)
    [a,b,c]=retta_per_2pti(xcbar(kk),ycbar(kk),x0,0);
    mOrto=-a/b/2;
    xmag=[xmag,cos(atan(mOrto))];
    ymag=[ymag,sin(atan(mOrto))];
    Br = [Br mat.LayerMag.Br(kk)];    % add one Br (or Air) block
end
temp.xc=xcbar;
temp.yc=ycbar;
temp.xmag=xmag;
temp.ymag=ymag;
temp.zmag=zeros(1,length(xcbar));

mat.LayerMag.Br = [mat.LayerMag.Br mat.LayerMag.Br];    % replicates Br for correct block assignation

% Salvataggio dei dati finali:
%%
hf=[r,B1k]-[B2k,Ar]; % rotor dx calculation in absolute term.
%%
geo.hf = hf;

temp.LowBarLength=LowBarLength;
temp.xc=xcbar;
temp.yc=ycbar;
% temp.xpont1=xpont1;
% temp.ypont1=ypont1;

temp.Bx0=Bx0;
temp.By0=By0;
temp.B1k=B1k;
temp.B2k=B2k;

temp.xpont=xpont;
temp.ypont=ypont;

% temp.RcirRib=rCirRib;
% geo.xpont2=xpont2;
% geo.ypont2=ypont2;

temp.xTraf1=xTraf1;
temp.xTraf2=xTraf2;
temp.yTraf1=yTraf1;
temp.yTraf2=yTraf2;
temp.arcLayer1=arcLayer1;
temp.arcLayer2=arcLayer2;
temp.arcLayTraf1=arcLayTraf1;
temp.arcLayTraf2=arcLayTraf2;
temp.XcRibTraf1=XcRibTraf1;
temp.YcRibTraf1=YcRibTraf1;
temp.XcRibTraf2=XcRibTraf2;
temp.YcRibTraf2=YcRibTraf2;
%% Punti per i ribs radiali
temp.XpontRadSx=XpontRadSx;
temp.YpontRadSx=YpontRadSx;
temp.XpontRadDx=XpontRadDx;
temp.YpontRadDx=YpontRadDx;
temp.XpontRadBarDx=XpontRadBarDx;
temp.YpontRadBarDx=YpontRadBarDx;
temp.XpontRadBarSx=XpontRadBarSx;
temp.YpontRadBarSx=YpontRadBarSx;

temp.xxD1k=xxD1k;
temp.yyD1k=yyD1k;
temp.xxD2k=xxD2k;
temp.yyD2k=yyD2k;

temp.XcBar1=XcBar1;
temp.YcBar1=YcBar1;
temp.XcBar2=XcBar2;
temp.YcBar2=YcBar2;
temp.XcBarLast_mean=XcBarLast_mean;
temp.YcBarLast_mean=YcBarLast_mean;

temp.xxB1k_mean=xxB1k_mean;
temp.yyB1k_mean=yyB1k_mean;
temp.xxB2k_mean=xxB2k_mean;
temp.yyB2k_mean=yyB2k_mean;
temp.xxB1k_mean2=xxB1k_mean2;
temp.yyB1k_mean2=yyB1k_mean2;
temp.xxB1k_mean3=xxB1k_mean3;
temp.yyB1k_mean3=yyB1k_mean3;
temp.arcLayer3=arcLayer3;
% temp.error_code=error_code;
% temp.r_fe=rfe;
% temp.x_fe=xFe;
% temp.y_fe=yFe;
temp.LastBarCurvatura=LastBarCurvatura;
geo.pontR=pont;


