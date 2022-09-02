
% Copyright 2018
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

function [temp,geo] = calc_ribs_rad_Seg(geo,mat,temp)

x0 = geo.x0;
r = geo.r;
l = geo.l;
p = geo.p;
pont0 = geo.pont0;              % Ponticelli al traferro (i ponticelli al traferro hanno lo spessore di un arco lungo pont0)
nlay = geo.nlay;
nmax = geo.nmax;                % Velocità max (rpm) per la valutazione della sollecitazione centrifuga più gravosa (-> ponticelli)
hc = geo.hc;
radial_ribs_split = geo.radial_ribs_split;  % flag to define if radial ribs is single or splitted
radial_ribs_eval  = geo.radial_ribs_eval;   % flag to select if ribs are automatically sized or not


pontRang       = geo.pontRang;
RotorFilletRad1 = geo.RotorFillet1;
RotorFilletRad2 = geo.RotorFillet2;

if ~strcmp(geo.RotType,'Seg')
    temp.flag_segV = zeros(1,nlay);
    temp.hcAngle = zeros(1,nlay);
end

B1k    = temp.B1k;
B2k    = temp.B2k;
flag_segV = temp.flag_segV;

rhoFE = mat.Rotor.kgm3;             % densità del ferro di rotore [kg/m3]
rhoPM = mat.LayerMag.kgm3;          % densità magneti [kg/m3]
sigma_max = mat.Rotor.sigma_max;    % snervamento materiale [MPa]

XpontRadBarSx = zeros(1,nlay);
YpontRadBarSx = zeros(1,nlay);
XpontRadBarDx = zeros(1,nlay);
YpontRadBarDx = zeros(1,nlay);
XpontRadDx    = zeros(1,nlay);
YpontRadDx    = zeros(1,nlay);
XpontRadSx    = zeros(1,nlay);
YpontRadSx    = zeros(1,nlay);
splitRib      = zeros(1,nlay);

XpontSplitBarSx = nan(2,geo.nlay);
YpontSplitBarSx = nan(2,geo.nlay);
XpontSplitSx    = nan(2,geo.nlay);
YpontSplitSx    = nan(2,geo.nlay);
XpontSplitDx    = nan(2,geo.nlay);
YpontSplitDx    = nan(2,geo.nlay);
XpontSplitBarDx = nan(2,geo.nlay);
YpontSplitBarDx = nan(2,geo.nlay);

xS01k = nan(1,nlay);
yS01k = nan(1,nlay);
xS02k = nan(1,nlay);
yS02k = nan(1,nlay);
xI01k = nan(1,nlay);
yI01k = nan(1,nlay);
xI02k = nan(1,nlay);
yI02k = nan(1,nlay);


%check if radial_ribs_split has just binary values
radial_ribs_split(radial_ribs_split>=0.5) = 1;
radial_ribs_split(radial_ribs_split<0.5) = 0;
radial_ribs_split(flag_segV==1) = 0;
geo.radial_ribs_split = radial_ribs_split;

RotorFilletRad2(RotorFilletRad2>hc/2)  = hc((RotorFilletRad2>hc/2))/2;
RotorFilletRad2(RotorFilletRad2<pont0) = pont0;
RotorFilletRad1(RotorFilletRad1>hc/2)  = hc((RotorFilletRad1>hc/2))/2;
RotorFilletRad1(RotorFilletRad1<pont0) = pont0;

% %limit Radial ribs Fillet
% RotorFilletRad1(RotorFilletRad1>hc/2)=round(hc((RotorFilletRad1>hc/2))/2,2);
% RotorFilletRad1(RotorFilletRad1<pont0)=pont0;

%check kOB with split ribs (angle feasible)
% pontRang(kOB<1 & pontRang==0)=-0.1;
% pontRang(pontRang>90/geo.p) = 90/geo.p;
% pontRang(pontRang<-85) = -85;

pontRang(pontRang>1) = 1;
pontRang(pontRang<0) = 0;
pontAng = pi-pontRang*pi/2*(p-1)/p;
pontAng(pontRang==0) = 0;
            

xpont  = temp.xpont;
ypont  = temp.ypont;
XpBar1 = temp.XpBar1;
YpBar1 = temp.YpBar1;
XpBar2 = temp.XpBar2;
YpBar2 = temp.YpBar2;
xxD1k  = temp.xxD1k;
yyD1k  = temp.yyD1k;
xxD2k  = temp.xxD2k;
yyD2k  = temp.yyD2k;

% % pont offset limit: down = half barrier / up = half barrier
% limDown = min([YpBar1;YpBar2])/2;
% xTmp1 = [xxD1k;xxD2k];
% yTmp1 = [yyD1k;yyD2k];
% xTmp2 = [XpBar1;YpBar1];
% yTmp2 = [YpBar1;YpBar2];
% limUp = min(((xTmp2-xTmp1)^2+(yTmp2-yTmp1)^2)^0.5);
% pontRoffset(pontRoffset>limDown) = limDown(pontRoffset>limDown);
% pontRoffset(pontRoffset>limUp) = limUp(pontRoffset>limUp);


% valid for SEG
% I-SEG: nearly correct
% Fluid: non correct

%%  Seg - calcolo area Ferro

% % check racc_pont: if the barriers are too small, racc_pont is set as hc/3
% racc_pont=racc_pont*ones(1,nlay);
% racc_pont(2*racc_pont>hc)=hc(2*racc_pont>hc)/3;

rTemp=abs((x0-xpont)+1j*ypont);
[xrot_traf,yrot_traf]=calc_intersezione_cerchi(r, rTemp, x0);
Dx=(r-xrot_traf(1))/5;  %per avere almeno 5 divisioni;
xcir_plot=[r:-Dx:r*cos(pi/2/p)];
ycir_plot=sqrt(r^2-xcir_plot.^2);
VectCir=find(xcir_plot>=xrot_traf(1));
x_ext_rot=xcir_plot(VectCir);
y_ext_rot=ycir_plot(VectCir);

A=[];
for ii=1:nlay
    if ii==1
        X=[B2k(ii), XpBar2(ii), xxD2k(ii),xpont(ii),fliplr(x_ext_rot)];
        Y=[0, YpBar2(ii), yyD2k(ii),ypont(ii),fliplr(y_ext_rot)];
        %         figure(100);hold on;fill(X,Y,'r');hold off;
        A(ii)=polyarea(X,Y);
        tmp_bary(ii,:)=centroid(X',Y');         %MarcoP
        clear X Y;
    else
        X=[B1k(ii-1), XpBar1(ii-1), xxD1k(ii-1),xpont(ii-1),xrot_traf(ii-1),xrot_traf(ii),xpont(ii),xxD2k(ii),XpBar2(ii),B2k(ii)];
        Y=[0, YpBar1(ii-1), yyD1k(ii-1),ypont(ii-1),yrot_traf(ii-1),yrot_traf(ii),ypont(ii),yyD2k(ii),YpBar2(ii),0];
        %         figure(100);hold on;fill(X,Y,'r'); hold off;
        A(ii)=polyarea(X,Y);
        tmp_bary(ii,:)=centroid(X',Y');         %MarcoP
        clear X Y;
    end
end
Afe = cumsum(A);

%MarcoP
rG(1)=tmp_bary(1,1);     % baricentro della regione di ferro sorretta dal ponticello a I
for ii=2:nlay
    rG(ii)=(Afe(ii-1)*rG(ii-1)+A(ii)*tmp_bary(ii,1))/Afe(ii);  % baricentri delle regioni di ferro sorrette dalle due U
end

M_Fe = 2*Afe*l * 1e-9 * rhoFE ;   % massa ferro appeso ai ponticelli

% Seg - calcolo area magnete
% barriera ad I
X_Ibarr=[xxD1k(1) xxD2k(1) B2k(1) B1k(1)];
Y_Ibarr=[yyD1k(1) yyD2k(1) 0 0];
areaI=2*polyarea(X_Ibarr,Y_Ibarr);  % area della barriera a I in cui inserire i magneti (escludo i due semicerchi all'estremità della barriera e gli eventuali ponticelli)
% % barriera ad U centrale
if nlay == 1
    area_barr_withoutPM = [areaI];
else
    % barriere ad U
    for kk=2:nlay
        areaUbars(kk-1) = (2*YpBar2(kk)*hc(kk))+(2*hc(kk)*abs((xxD2k(kk)+i*yyD2k(kk))-(XpBar2(kk)+i*YpBar2(kk)))); %area della barriera a U in cui inserire i magneti (escludo i due semicerchi all'estremità della barriera)
    end
    area_barr_withoutPM = [areaI areaUbars];
end
rG_PM = (B1k+B2k)/2;    % vettore con i baricentri dei magneti
area_barr_withPM=area_barr_withoutPM;
A_PM=cumsum(area_barr_withPM);
M_PM=A_PM*rhoPM*1e-9*2*l;

% Seg - calcolo e disegno ponticelli
F_centrifuga = (M_Fe+M_PM) .* rG/1000 *  (nmax * pi/30)^2;

if radial_ribs_eval == 0
    pont = F_centrifuga/(sigma_max * l);    % mm
else
    pont = geo.pontR;
end

if ~flag_segV
    pont(pont>(YpBar2-pont0)) = YpBar2(pont>(YpBar2-pont0))-pont0*2;
else
    pont(pont>(yyD2k-pont0)) = yyD2k(pont>(yyD2k-pont0))-pont0*2;
end

% do not use split ribs if barrier is small
radial_ribs_split(YpBar2<pont*1.5) = 0;


%% Ribs definition

flagPosSplitDx  = nan(2,nlay);
flagPosSplitSx  = nan(2,nlay);
% le due variabili di sopra definiscono dove iniziano i ponticelli e vanno
% come le variabili splitRibs. Se 0, sono nella parte centrale, se 1 sono
% nella parte esterna

% Prima controllo lo spessore minimo dei ponticelli, anche in relazione a
% split ribs o no.
% NB: posso decidere qui se eliminare i ribs splittati se ponticello troppo
% piccolo
for ii=1:nlay
    if radial_ribs_split(ii)==0
        if pont(ii)<pont0
            pont(ii)=0;
        end
    else
        if pont(ii)<pont0
            pont(ii)=0;
        elseif pont(ii)<2*pont0
            pont(ii) = 2*pont0;
        end
    end
end

% Disegno ponticelli con raccordi

for ii=1:nlay
    if pont(ii)==0
        XpontRadBarSx(ii) = B1k(ii);
        YpontRadBarSx(ii) = 0;
        XpontRadBarDx(ii) = B2k(ii);
        YpontRadBarDx(ii) = 0;
        XpontRadDx(ii)    = NaN;
        YpontRadDx(ii)    = 0;
        XpontRadSx(ii)    = NaN;
        YpontRadSx(ii)    = 0;

        XpontSplitBarSx(1,ii) = XpBar1(ii);
        XpontSplitBarSx(2,ii) = XpBar1(ii);
        YpontSplitBarSx(1,ii) = YpBar1(ii);
        YpontSplitBarSx(2,ii) = YpBar1(ii);
        XpontSplitBarDx(1,ii) = XpBar2(ii);
        XpontSplitBarDx(2,ii) = XpBar2(ii);
        YpontSplitBarDx(1,ii) = YpBar2(ii);
        YpontSplitBarDx(2,ii) = YpBar2(ii);
    else
        if radial_ribs_split(ii)==0
            if flag_segV(ii)==1
                a0 = 0;
                b0 = 1;
                c0 = -pont(ii)/2;
                % Raccordo lato sx (lato 1)
                [a1,b1,c1] = retta_per_2pti(B1k(ii),0,xxD1k(ii),yyD1k(ii));
                % Calcolo le rette offset del raggio di raccordo
                [a0o,b0o,c0o] = calc_retta_offset(a0,b0,c0,-RotorFilletRad1(ii));
                [a1o,b1o,c1o] = calc_retta_offset(a1,b1,c1,+RotorFilletRad1(ii));
                [xS01k(ii),yS01k(ii)] = intersezione_tra_rette(a0o,b0o,c0o,a1o,b1o,c1o);
                % Trovo le proiezioni del centro del raccordo sulle 2 rette, per definire dove inizia e finisce l'arco
                [XpontRadBarSx(ii),YpontRadBarSx(ii)]       = proiezione_punto_retta(a1,b1,c1,xS01k(ii),yS01k(ii));
                [XpontRadSx(ii),YpontRadSx(ii)] = proiezione_punto_retta(a0,b0,c0,xS01k(ii),yS01k(ii));
                % Raccordo lato dx (lato 2)
                [a2,b2,c2] = retta_per_2pti(B2k(ii),0,xxD2k(ii),yyD2k(ii));
                % Calcolo le rette offset del raggio di raccordo
                [a0o,b0o,c0o] = calc_retta_offset(a0,b0,c0,-RotorFilletRad2(ii));
                [a2o,b2o,c2o] = calc_retta_offset(a2,b2,c2,-RotorFilletRad2(ii));
                [xS02k(ii),yS02k(ii)] = intersezione_tra_rette(a0o,b0o,c0o,a2o,b2o,c2o);
                % Trovo le proiezioni del centro del raccordo sulle 2 rette, per definire dove inizia e finisce l'arco
                [XpontRadBarDx(ii),YpontRadBarDx(ii)]       = proiezione_punto_retta(a2,b2,c2,xS02k(ii),yS02k(ii));
                [XpontRadDx(ii),YpontRadDx(ii)] = proiezione_punto_retta(a0,b0,c0,xS02k(ii),yS02k(ii));
                XpBar1(ii) = XpontRadBarSx(ii);
                YpBar1(ii) = YpontRadBarSx(ii);
                XpBar2(ii) = XpontRadBarDx(ii);
                YpBar2(ii) = YpontRadBarDx(ii);

                XpontSplitBarSx(1,ii) = XpBar1(ii);
                XpontSplitBarSx(2,ii) = XpBar1(ii);
                YpontSplitBarSx(1,ii) = YpBar1(ii);
                YpontSplitBarSx(2,ii) = YpBar1(ii);
                XpontSplitBarDx(1,ii) = XpBar2(ii);
                XpontSplitBarDx(2,ii) = XpBar2(ii);
                YpontSplitBarDx(1,ii) = YpBar2(ii);
                YpontSplitBarDx(2,ii) = YpBar2(ii);

            else
                XpontRadSx(ii)    = B1k(ii)+RotorFilletRad1(ii);
                YpontRadSx(ii)    = pont(ii)/2;
                XpontRadDx(ii)    = B2k(ii)-RotorFilletRad2(ii);
                YpontRadDx(ii)    = pont(ii)/2;
                XpontRadBarSx(ii) = B1k(ii);
                YpontRadBarSx(ii) = pont(ii)/2+RotorFilletRad1(ii);
                XpontRadBarDx(ii) = B2k(ii);
                YpontRadBarDx(ii) = pont(ii)/2+RotorFilletRad2(ii);
                xS01k(ii)         = XpontRadSx(ii);
                yS01k(ii)         = YpontRadBarSx(ii);
                xS02k(ii)         = XpontRadDx(ii);
                yS02k(ii)         = YpontRadBarDx(ii);

                XpontSplitBarSx(1,ii) = XpBar1(ii);
                XpontSplitBarSx(2,ii) = XpBar1(ii);
                YpontSplitBarSx(1,ii) = YpBar1(ii);
                YpontSplitBarSx(2,ii) = YpBar1(ii);
                XpontSplitBarDx(1,ii) = XpBar2(ii);
                XpontSplitBarDx(2,ii) = XpBar2(ii);
                YpontSplitBarDx(1,ii) = YpBar2(ii);
                YpontSplitBarDx(2,ii) = YpBar2(ii);


            end
        else
            % Nomenclatura:
            % - retta p = asse ponticello
            % - retta 0 = lato ponticello
            % - rette 1,2 = lati barriera
            % - rette _o = offset (per calcolo raccordo)
            ap = tan(pontAng(ii));
            bp = -1;
            cp = (YpBar2(ii))-ap*XpBar2(ii);
            % parte inferiore del ponticello (seconda riga variabile)
            [a0,b0,c0] = calc_retta_offset(ap,bp,cp,-pont(ii)/4);
            % lato dx
            [a2,b2,c2] = retta_per_2pti(B2k(ii),0,XpBar2(ii),YpBar2(ii));
            flagPosSplitDx(2,ii) = 0;
            % calcolo raccordo: rette offset per centro e proiezioni del centro
            % sulle rette da raccordare
            [a0o,b0o,c0o] = calc_retta_offset(a0,b0,c0,-RotorFilletRad2(ii));
            [a2o,b2o,c2o] = calc_retta_offset(a2,b2,c2,-RotorFilletRad2(ii));
            [xI02k(ii),yI02k(ii)] = intersezione_tra_rette(a0o,b0o,c0o,a2o,b2o,c2o);
            [XpontSplitBarDx(2,ii),YpontSplitBarDx(2,ii)] = proiezione_punto_retta(a2,b2,c2,xI02k(ii),yI02k(ii));
            [XpontSplitDx(2,ii),YpontSplitDx(2,ii)] = proiezione_punto_retta(a0,b0,c0,xI02k(ii),yI02k(ii));
            % lato sx
            [a0o,b0o,c0o] = calc_retta_offset(a0,b0,c0,-RotorFilletRad1(ii));
            if YpBar1(ii)<(-B1k(ii)*a0o-c0o)/b0o
                [a1,b1,c1] = retta_per_2pti(xxD1k(ii),yyD1k(ii),XpBar1(ii),YpBar1(ii));
                flagPosSplitSx(2,ii) = 1;
                [a1o,b1o,c1o] = calc_retta_offset(a1,b1,c1,-RotorFilletRad1(ii));
            else
                [a1,b1,c1] = retta_per_2pti(B1k(ii),0,XpBar1(ii),YpBar1(ii));
                flagPosSplitSx(2,ii) = 0;
                [a1o,b1o,c1o] = calc_retta_offset(a1,b1,c1,+RotorFilletRad1(ii));
            end
            [xI01k(ii),yI01k(ii)] = intersezione_tra_rette(a0o,b0o,c0o,a1o,b1o,c1o);
            [XpontSplitBarSx(2,ii),YpontSplitBarSx(2,ii)] = proiezione_punto_retta(a1,b1,c1,xI01k(ii),yI01k(ii));
            [XpontSplitSx(2,ii),YpontSplitSx(2,ii)] = proiezione_punto_retta(a0,b0,c0,xI01k(ii),yI01k(ii));
            % lato superiore del ponticello (prima riga della variabile)
            [a0,b0,c0] = calc_retta_offset(ap,bp,cp,+pont(ii)/4);
            % raccordo lato dx
            [a2,b2,c2] = retta_per_2pti(xxD2k(ii),yyD2k(ii),XpBar2(ii),YpBar2(ii));
            flagPosSplitDx(1,ii) = 1;
            % calcolo rette offset per centro raccordo
            [a0o,b0o,c0o] = calc_retta_offset(a0,b0,c0,+RotorFilletRad2(ii));
            [a2o,b2o,c2o] = calc_retta_offset(a2,b2,c2,+RotorFilletRad2(ii));
            if xxD2k(ii)==XpBar2(ii) % per prima barriera (se cicciona), non ho il tratto rettilineo --> no raccordo
                xS02k(ii)             = xxD2k(ii);
                yS02k(ii)             = yyD2k(ii);
                XpontSplitDx(1,ii)    = xxD2k(ii);
                YpontSplitDx(1,ii)    = yyD2k(ii);
                XpontSplitBarDx(1,ii) = xxD2k(ii);
                YpontSplitBarDx(1,ii) = yyD2k(ii);
            else
                [xS02k(ii),yS02k(ii)] = intersezione_tra_rette(a0o,b0o,c0o,a2o,b2o,c2o);
                [XpontSplitBarDx(1,ii),YpontSplitBarDx(1,ii)] = proiezione_punto_retta(a2,b2,c2,xS02k(ii),yS02k(ii));
                [XpontSplitDx(1,ii),YpontSplitDx(1,ii)] = proiezione_punto_retta(a0,b0,c0,xS02k(ii),yS02k(ii));
            end
            % raccordo lato sx
            [a0o,b0o,c0o] = calc_retta_offset(a0,b0,c0,+RotorFilletRad1(ii));
            if YpBar1(ii)<(-B1k(ii)*a0o-c0o)/b0o
                [a1,b1,c1] = retta_per_2pti(xxD1k(ii),yyD1k(ii),XpBar1(ii),YpBar1(ii));
                flagPosSplitSx(1,ii) = 1;
                [a1o,b1o,c1o] = calc_retta_offset(a1,b1,c1,-RotorFilletRad1(ii));
            else
                [a1,b1,c1] = retta_per_2pti(B1k(ii),0,XpBar1(ii),YpBar1(ii));
                flagPosSplitSx(1,ii) = 0;
                [a1o,b1o,c1o] = calc_retta_offset(a1,b1,c1,+RotorFilletRad1(ii));
            end
            % calcolo offset e raccordo
            [xS01k(ii),yS01k(ii)] = intersezione_tra_rette(a0o,b0o,c0o,a1o,b1o,c1o);
            [XpontSplitBarSx(1,ii),YpontSplitBarSx(1,ii)] = proiezione_punto_retta(a1,b1,c1,xS01k(ii),yS01k(ii));
            [XpontSplitSx(1,ii),YpontSplitSx(1,ii)] = proiezione_punto_retta(a0,b0,c0,xS01k(ii),yS01k(ii));

            % check if fillet is small enough
            if YpontSplitBarDx(1,ii)>yyD2k(ii)
                XpontSplitDx(1,ii)    = xxD2k(ii);
                YpontSplitDx(1,ii)    = yyD2k(ii);
                XpontSplitBarDx(1,ii) = xxD2k(ii);
                YpontSplitBarDx(1,ii) = xxD2k(ii);
                xS02k(ii)             = xxD2k(ii);
                yS02k(ii)             = yyD2k(ii);
            end
            XpontRadBarDx(ii) = B2k(ii);
            YpontRadBarDx(ii) = 0;
            XpontRadDx(ii)    = NaN;
            YpontRadDx(ii)    = NaN;
            XpontRadBarSx(ii) = B1k(ii);
            YpontRadBarSx(ii) = 0;
            XpontRadSx(ii)    = NaN;
            YpontRadSx(ii)    = NaN;
        end
    end
end


% Assign output values to geo and temp
temp.XpontRadDx    = XpontRadDx;
temp.YpontRadDx    = YpontRadDx;
temp.XpontRadSx    = XpontRadSx;
temp.YpontRadSx    = YpontRadSx;
temp.XpontRadBarDx = XpontRadBarDx;
temp.XpontRadBarSx = XpontRadBarSx;
temp.YpontRadBarDx = YpontRadBarDx;
temp.YpontRadBarSx = YpontRadBarSx;
temp.splitRib      = splitRib;

temp.XpontSplitBarSx = XpontSplitBarSx;
temp.YpontSplitBarSx = YpontSplitBarSx;
temp.XpontSplitBarDx = XpontSplitBarDx;
temp.YpontSplitBarDx = YpontSplitBarDx;
temp.XpontSplitDx    = XpontSplitDx;
temp.YpontSplitDx    = YpontSplitDx;
temp.XpontSplitSx    = XpontSplitSx;
temp.YpontSplitSx    = YpontSplitSx;

temp.xS01k = xS01k;
temp.yS01k = yS01k;
temp.xS02k = xS02k;
temp.yS02k = yS02k;

temp.xI01k = xI01k;
temp.yI01k = yI01k;
temp.xI02k = xI02k;
temp.yI02k = yI02k;

temp.XpBar1 = XpBar1;
temp.YpBar1 = YpBar1;
temp.XpBar2 = XpBar2;
temp.YpBar2 = YpBar2;
temp.xxD1k  = xxD1k;
temp.yyD1k  = yyD1k;
temp.xxD2k  = xxD2k;
temp.yyD2k  = yyD2k;



% temp.flag_cent    = flag_cent;
% temp.flag_ext     = flag_ext;
% temp.flag_shift   = flag_shift;
% temp.flag_shiftUP = flag_shiftUP;

temp.flagPosSplitDx = flagPosSplitDx;
temp.flagPosSplitSx = flagPosSplitSx;

geo.pontRang = pontRang;
% geo.pontRoffset = pontRoffset;
geo.pontR = pont;
geo.RotorFillet1 = RotorFilletRad1;
geo.RotorFillet2 = RotorFilletRad2;

geo.radial_ribs_split = radial_ribs_split;





