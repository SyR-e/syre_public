
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

function [temp,geo] = calc_ribs_rad_fun(geo,mat,temp)

% NEW: function version of script calc_ribs_rad
% OK for CIRCULAR and SEG
% Other geometries still use the old script
% after fixing alla geometries, the function name must return to "calc_ribs_rad"

x0 = geo.x0;
r = geo.r;
l = geo.l;
p = geo.p;
pont0 = geo.pont0;              % Ponticelli al traferro (i ponticelli al traferro hanno lo spessore di un arco lungo pont0)
racc_pont = geo.racc_pont;      % racc_pont=1*pont0 <- per i ponticelli radiali.
nlay = geo.nlay;
nmax = geo.nmax;                % Velocità max (rpm) per la valutazione della sollecitazione centrifuga più gravosa (-> ponticelli)
hc = geo.hc;
radial_ribs_split = geo.radial_ribs_split;  % flag to define if radial ribs is single or splitted
radial_ribs_eval  = geo.radial_ribs_eval;   % flag to select if ribs are automatically sized or not

pontRang       = geo.pontRang;
pontRoffset    = geo.pontRoffset;
RotorFilletRad = geo.RotorFillet1;
RotorFilletRad1 = geo.RotorFillet1;
RotorFilletRad2 = geo.RotorFillet2;
kOB             = geo.kOB;
delta_FBS       = geo.delta_FBS;

if ~strcmp(geo.RotType,'Seg')
    temp.flag_segV = zeros(1,nlay);
    temp.hcAngle = zeros(1,nlay);
end

B1k    = temp.B1k;
B2k    = temp.B2k;
flag_segV = temp.flag_segV;
hcAngle   = temp.hcAngle;

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

xS01k           = zeros(1,nlay);
yS01k           = zeros(1,nlay);
xS02k           = zeros(1,nlay);
yS02k           = zeros(1,nlay);

xI01k = zeros(1,nlay);
yI01k = zeros(1,nlay);
xI02k = zeros(1,nlay);
yI02k = zeros(1,nlay);

flag_cent    = zeros(1,nlay);
flag_ext     = zeros(1,nlay);
flag_shift   = zeros(1,nlay);
flag_shiftUP = zeros(1,nlay);


%check if radial_ribs_split has just binary values
radial_ribs_split(radial_ribs_split>=0.5) = 1;
radial_ribs_split(radial_ribs_split<0.5) = 0;
radial_ribs_split(flag_segV==1) = 0;
geo.radial_ribs_split = radial_ribs_split;

RotorFilletRad2(RotorFilletRad2>hc/2)=round(hc((RotorFilletRad2>hc/2))/2,2);
RotorFilletRad1(RotorFilletRad1<pont0)=pont0;RotorFilletRad1(RotorFilletRad1>hc/2)=round(hc((RotorFilletRad1>hc/2))/2,2);
RotorFilletRad1(RotorFilletRad1<pont0)=pont0;

%limit Radial ribs Fillet
RotorFilletRad(RotorFilletRad>hc/2)=round(hc((RotorFilletRad>hc/2))/2,2);
RotorFilletRad(RotorFilletRad<pont0)=pont0;

%check kOB with split ribs (angle feasible)
pontRang(kOB<1 & pontRang==0)=-0.1;

switch geo.RotType
    case 'Circular'
        
        % Valutazione ponticelli radiali:
        r_all = [];                                                                 % Nel vettore r_all si memorizzano invece le distanze, sempre rispetto al
        % Starting and ending point of the flux barrier...
        for jj = 1:nlay                                                             % che individuano l'inizio e la fine delle nlay barriere.
            r_all = [r_all x0-B2k(jj) x0-B1k(jj)];
        end
        % median point of flux barrier...
        rbeta=(r_all(1:2:end)+r_all(2:2:end))./2;
        % Posizione banane su asse q (convenzione assi VAGATI)
        XBanqdx=x0-r_all(1:2:end);
        XBanqsx=x0-r_all(2:2:end);
        
        A = zeros(1,nlay);  % carrier Area (Fe)
        Ab = A;             % barrier Area (Air or PM)
        tmp_bary = A;       % carrier center of gravity
        tmp_bary_b = A;     % barrier center of gravity
        rG = A;             % aggregate center of gravity (nlay equivalent masses, Fe + PM)
        
        % carriers (Fe): Area and Center of Gravity
        for j = 1 : nlay
            if j == 1
                % first carrier
                [x,y] = calc_intersezione_cerchi(r,r_all(1),x0);
                dist = x0 - x;
                theta = atan2(y,dist);
                theta2 = atan2(y,x);
                A1 = 0.5*(r^2*theta2 - abs(0.5*det([x y 1; 0 0 1; x -y 1])));
                A2 = 0.5*(r_all(j)^2*theta - abs(0.5*det([x y 1; x0 0 1; x -y 1])));
                A(j) = A1+A2;
                bary1 = (2*r*(sin(theta2))^3)/(3*(theta2-sin(theta2)*cos(theta2)));
                bary2 = x0-(2*r_all(1)*(sin(theta))^3)/(3*(theta-sin(theta)*cos(theta)));
                tmp_bary(j) = (2*A1*bary1 + 2*A2*bary2)/(2*(A1+A2));
            else
                [x,y] = calc_intersezione_cerchi(r,r_all(2*j-1),x0);
                dist = x0 - x;
                theta = atan2(y,dist);
                A(j) = (r_all(2*j-1)^2 - r_all(2*j-2)^2)*theta*0.5;
                tmp_bary(j) = x0 - (2/3*sin(theta)/theta*(r_all(2*j-1)^3 - r_all(2*j-2)^3)/(r_all(2*j-1)^2 - r_all(2*j-2)^2));
            end
        end
        
        % barriers (Air or PM): Area and Center of Gravity
        for j = 1 : nlay
            [x,y] = calc_intersezione_cerchi(r,r_all(2*j),x0);
            dist = x0 - x;
            theta = atan2(y,dist);
            Ab(j) = (r_all(2*j)^2 - r_all(2*j-1)^2)*theta*0.5;
            tmp_bary_b(j) = x0 - (2/3*sin(theta)/theta*(r_all(2*j)^3 - r_all(2*j-1)^3)/(r_all(2*j)^2 - r_all(2*j-1)^2));
        end
        
        % islands as seen by the bridges
        Afe = cumsum(A);
        Abarr = cumsum(Ab);
        
        % mass and center of gravity of interest for bridges
        if (mean(unique(mat.LayerMag.Br)) > 0)
            % Fe + PM
            mass = Afe*rhoFE+Abarr*rhoPM; %*geo.BarFillFac;
            %rG(1) = (tmp_bary(1)*Afe(1)*rhoFE+tmp_bary_b(1)*Abarr(1)*rhoPM*geo.BarFillFac)/(mass(1));
            rG(1) = (tmp_bary(1)*Afe(1)*rhoFE+tmp_bary_b(1)*Abarr(1)*rhoPM)/(mass(1));
            for ii = 2 : nlay
                %rG(ii)=(mass(ii-1)*rG(ii-1)+tmp_bary(ii)*A(ii)*rhoFE+tmp_bary_b(ii)*Ab(ii)*rhoPM*geo.BarFillFac)/mass(ii);  % baricentri delle regioni di ferro sorrette dalle due U
                rG(ii)=(mass(ii-1)*rG(ii-1)+tmp_bary(ii)*A(ii)*rhoFE+tmp_bary_b(ii)*Ab(ii)*rhoPM)/mass(ii);  % baricentri delle regioni di ferro sorrette dalle due U
            end
        else
            % Fe + Air
            mass = Afe*rhoFE;
            rG(1) = tmp_bary(1);
            for ii = 2 : nlay
                rG(ii)=(Afe(ii-1)*rG(ii-1)+A(ii)*tmp_bary(ii))/Afe(ii);  % baricentri delle regioni di ferro sorrette dalle due U
            end
        end
        
        M_tot = mass * 2 * l * 1e-9;   % kg
        F_centrifuga = M_tot .* rG/1000 * (nmax * pi/30)^2;
        
        if radial_ribs_eval == 0
            pont = F_centrifuga/(sigma_max * l);    % mm
        else
            pont = geo.pontR;                        % input from GUI
        end
        
        pont(pont<pont0) = 0;   % pont dimension >= pont0
        
        racc_pont = racc_pont*ones(1,nlay);
        racc_pont(2*racc_pont>hc) = hc(2*racc_pont>hc)/3;
        
        %%% Radial Ribs Node Coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for jj = 1 : 2 : length(r_all)
            x10 = x0 - r_all(jj) - racc_pont(ceil(jj/2));   % (Ponticello radiale_vertice in alto a dx)=(Punto in basso a dx del raccordo))
            x11 = x0 - r_all(jj+1) + racc_pont(ceil(jj/2)); % (Ponticello radiale_vertice in alto a sx)=(Punto in basso a sx del raccordo))
            if (x11 < x10) && (pont(ceil(jj/2))>0)
                % Se i due punti non si incrociano (raccordo non troppo
                % grande rispetto alla larghezza barriera)
                hpont = pont(ceil(jj/2));
                y10 = hpont/2;
                YpontRadBarDx(ceil(jj/2)) = y10 + racc_pont(ceil(jj/2));
                [x_temp,~] = intersezione_retta_circonferenza(x0,0,r_all(jj),0,YpontRadBarDx(ceil(jj/2)));
                XpontRadBarDx(ceil(jj/2)) = 2*x0 - x_temp;
                YpontRadBarSx(ceil(jj/2)) = y10 + racc_pont(ceil(jj/2));
                [x_temp,~] = intersezione_retta_circonferenza(x0,0,r_all(jj+1),0,YpontRadBarSx(ceil(jj/2)));
                XpontRadBarSx(ceil(jj/2)) = 2*x0 - x_temp;
                XpontRadDx(ceil(jj/2))  = XpontRadBarDx(ceil(jj/2))  - racc_pont(ceil(jj/2));
                YpontRadDx(ceil(jj/2))  = y10;
                XpontRadSx(ceil(jj/2))  = XpontRadBarSx(ceil(jj/2))  + racc_pont(ceil(jj/2));
                YpontRadSx(ceil(jj/2))  = y10;
                
            else
                % Se invece i punti 10 e 11 si incrociano (ovvero raccordo troppo grande rispetto)
                % alla larghezza barriera), non faccio il ponticello e disegno solo una linea
                hpont = 0;
                XpontRadBarSx(ceil(jj/2)) = XBanqsx(ceil(jj/2));
                YpontRadBarSx(ceil(jj/2)) = 0;
                XpontRadBarDx(ceil(jj/2)) = XBanqdx(ceil(jj/2));
                YpontRadBarDx(ceil(jj/2)) = 0;
                XpontRadDx(ceil(jj/2)) = NaN;
                YpontRadDx(ceil(jj/2)) = 0;
                XpontRadSx(ceil(jj/2)) = NaN;
                YpontRadSx(ceil(jj/2)) = 0;
            end
        end
        
        geo.pontR = pont;
        geo.r_all = r_all;
        
        % Points for radial ribs
        temp.XpontRadDx    = XpontRadDx;
        temp.YpontRadDx    = YpontRadDx;
        temp.XpontRadSx    = XpontRadSx;
        temp.YpontRadSx    = YpontRadSx;
        temp.XpontRadBarDx = XpontRadBarDx;
        temp.XpontRadBarSx = XpontRadBarSx;
        temp.YpontRadBarDx = YpontRadBarDx;
        temp.YpontRadBarSx = YpontRadBarSx;
        temp.splitRib      = splitRib;
        
        temp.XpontSplitBarSx = XpontSplitBarSx; % coordinates of split inner ribs
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
        
        temp.flag_cent    = flag_cent;
        temp.flag_ext     = flag_ext;
        temp.flag_shift   = flag_shift;
        temp.flag_shiftUP = flag_shiftUP;
        
        geo.pontRang=pontRang;
        geo.pontRoffset=pontRoffset;
        geo.RotorFillet(2,:)=RotorFilletRad;
        geo.pontR=pont;
        
    otherwise
        
        % valid for SEG
        % I-SEG: nearly correct
        % Fluid: non correct
        
        %  Seg - calcolo area Ferro
        
        
        % check racc_pont: if the barriers are too small, racc_pont is set as hc/3
        racc_pont=racc_pont*ones(1,nlay);
        racc_pont(2*racc_pont>hc)=hc(2*racc_pont>hc)/3;
        
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
        
        %         XpontSplitBarSx = nan(2,nlay);
        %         YpontSplitBarSx = nan(2,nlay);
        %         XpontSplitBarDx = nan(2,nlay);
        %         YpontSplitBarDx = nan(2,nlay);
        %         XpontSplitDx    = nan(2,nlay);
        %         YpontSplitDx    = nan(2,nlay);
        %         XpontSplitSx    = nan(2,nlay);
        %         YpontSplitSx    = nan(2,nlay);
        
        
        
        %         switch radial_ribs_split
        %
        %             case 0 % single rib
        
        if any(radial_ribs_split==0)
            for jj=1:nlay
                if (pont(jj) < pont0) % non disegno i ponticelli radiali il cui spessore è minore della tolleranza di lavorazione per gli altri tipi di rotore
                    pont(jj)=0;
                end
            end
            
            hpont=pont/2;
            %rac_pont=abs(B1k-B2k)/4;
            %rac_pont=pont0;
            
            if any(flag_segV==1)
                ii = (flag_segV==1);
                
                %                x_delta(ii) = pont(ii)/2./tan(pi/2/p+delta_FBS/2+hcAngle(ii));
                %                YpontRadSx(ii) = pont(ii)/2;
                %                YpontRadDx(ii) = pont(ii)/2;
                %                XpontRadSx(ii) = B1k(ii)+x_delta(ii)+RotorFilletRad1(ii);
                %                XpontRadDx(ii) = B2k(ii)+x_delta(ii)-RotorFilletRad2(ii);
                %
                %                x1(ii) = B1k(ii)+x_delta(ii);
                %                x2(ii) = B2k(ii)+x_delta(ii);
                %
                %
                %                YpontRadBarSx(ii) = pont(ii)/2+ RotorFilletRad1(ii).*cos(pi/2/p+delta_FBS/2+hcAngle(ii));
                %                YpontRadBarDx(ii) = pont(ii)/2+ RotorFilletRad2(ii).*cos(pi/2/p+delta_FBS/2+hcAngle(ii));
                %                XpontRadBarSx(ii) = B1k(ii)+ RotorFilletRad1(ii).*sin(pi/2/p+delta_FBS/2+hcAngle(ii));
                %                XpontRadBarDx(ii) = B2k(ii)+ RotorFilletRad2(ii).*sin(pi/2/p+delta_FBS/2+hcAngle(ii));
                %
                %                [a b c]=retta_per_2pti(B1k,0,XpontRadBarSx,YpontRadBarSx);
                %                [m q] = retta_abc2mq(a, b, c);
                %
                %                q = YpontRadSx - m.*XpontRadSx;
                %
                %                yS01k(ii) = YpontRadBarSx(ii);
                %                xS01k(ii) = (yS01k(ii)-q(ii))./m(ii);
                %
                %                [a b c]=retta_per_2pti(B2k,0,XpontRadBarDx,YpontRadBarDx);
                %                [m q] = retta_abc2mq(a, b, c);
                %
                %                q = YpontRadDx - m.*XpontRadDx;
                %
                %                yS02k(ii) = YpontRadBarDx(ii);
                %                xS02k(ii) = (yS02k(ii)-q(ii))./m(ii);
                
                a0 = 0;
                b0 = 1;
                c0 = -pont(ii)/2;
                [a1,b1,c1] = retta_per_2pti(B1k(ii),0,xxD1k(ii),yyD1k(ii));
                [a2,b2,c2] = retta_per_2pti(B2k(ii),0,xxD2k(ii),yyD2k(ii));
                
                [m0,q0,~]=retta_abc2mq(a0,b0,c0);
                [m1,q1,~]=retta_abc2mq(a1,b1,c1);
                [m2,q2,~]=retta_abc2mq(a2,b2,c2);
                
                yp1 = pont(ii)/2;
                xp1 = (yp1-q1)./m1;
                
                yp2 = yp1;
                xp2 = (yp2-q2)./m2;
                
                XpontRadSx(ii) = xp1+RotorFilletRad1(ii);
                YpontRadSx(ii) = yp1;
                
                XpontRadDx(ii) = xp2-RotorFilletRad2(ii);
                YpontRadDx(ii) = yp2;
                
                a = (1+m1.*m1);
                b = -2*xp1+2*m1.*(q1-yp1);
                c = -RotorFilletRad1(ii).*RotorFilletRad1(ii)+xp1.^2+(q1-yp1).^2;
                XpontRadBarSx(ii) = (-b+sqrt(b.^2-4*a.*c))./(2*a);
                YpontRadBarSx(ii) = m1.*XpontRadBarSx(1,(ii))+q1;
                
                a = (1+m2.*m2);
                b = -2*xp2+2*m2.*(q2-yp2);
                c = -RotorFilletRad2(ii).*RotorFilletRad2(ii)+xp2.^2+(q2-yp2).^2;
                XpontRadBarDx(ii) = (-b+sqrt(b.^2-4*a.*c))./(2*a);
                YpontRadBarDx(ii) = m2.*XpontRadBarDx(1,(ii))+q2;
                
                m1perp = -(m1.^-1);
                q1perp = YpontRadBarSx(ii)-m1perp.*XpontRadBarSx(ii);
                
                m2perp = -(m2.^-1);
                q2perp = YpontRadBarDx(ii)-m2perp.*XpontRadBarDx(ii);
                
                xS01k(ii) = XpontRadSx(ii);
                yS01k(ii) = m1perp.*xS01k(ii)+q1perp;
                
                xS02k(ii) = XpontRadDx(ii);
                yS02k(ii) = m2perp.*xS02k(ii)+q2perp;
                
                XpBar1(ii) = XpontRadBarSx(ii);
                YpBar1(ii) = YpontRadBarSx(ii);
                
                XpBar2(ii) = XpontRadBarDx(ii);
                YpBar2(ii) = YpontRadBarDx(ii);
                
         
                
                
            end
            
            %for ii=1:nlay
            if any(any(radial_ribs_split==0) & hpont>0 & flag_segV ==0)
                %if hpont(ii)>0
                ii=(radial_ribs_split==0 & hpont>0 & flag_segV ==0);
                
                XpontRadSx(ii)=B1k(ii)+RotorFilletRad1(ii);
                YpontRadSx(ii)=pont(ii)/2;
                
                XpontRadDx(ii)=B2k(ii)-RotorFilletRad2(ii);
                YpontRadDx(ii)=pont(ii)/2;
                
                XpontRadBarSx(ii)=B1k(ii);
                YpontRadBarSx(ii)=pont(ii)/2+RotorFilletRad1(ii);
                
                XpontRadBarDx(ii)=B2k(ii);
                YpontRadBarDx(ii)=pont(ii)/2+RotorFilletRad2(ii);
                
                xS01k(ii)=XpontRadSx(ii);
                yS01k(ii)=YpontRadBarSx(ii);
                
                xS02k(ii)=XpontRadDx(ii);
                yS02k(ii)=YpontRadBarDx(ii);
            end
            
            if any(any(radial_ribs_split==0) & hpont<=0 & flag_segV ==0)
                ii=(radial_ribs_split==0 & hpont<=0 & flag_segV ==0);
                XpontRadBarSx(ii)=B1k(ii);
                YpontRadBarSx(ii)=0;
                XpontRadBarDx(ii)=B2k(ii);
                YpontRadBarDx(ii)=0;
                XpontRadDx(ii)=NaN;
                YpontRadDx(ii)=0;
                XpontRadSx(ii)=NaN;
                YpontRadSx(ii)=0;
            end
            
            ii=radial_ribs_split==0;
            
            XpontSplitBarSx(1,ii)=B1k(ii);
            XpontSplitBarSx(2,ii)=B1k(ii);
            YpontSplitBarSx(1,ii)=YpBar1(ii);
            YpontSplitBarSx(2,ii)=YpontSplitBarSx(1,ii);
            XpontSplitBarDx(1,ii)=B2k(ii);
            XpontSplitBarDx(2,ii)=B2k(ii);
            YpontSplitBarDx(1,ii)=YpBar2(ii);
            YpontSplitBarDx(2,ii)=YpontSplitBarDx(1,ii);
            
            jj = (ii & flag_segV);
            XpontSplitBarSx(1,jj)=XpBar1(jj);
            XpontSplitBarSx(2,jj)=XpBar1(jj);
            YpontSplitBarSx(1,jj)=YpBar1(jj);
            YpontSplitBarSx(2,jj)=YpBar1(jj);
            XpontSplitBarDx(1,jj)=XpBar2(jj);
            XpontSplitBarDx(2,jj)=XpBar2(jj);
            YpontSplitBarDx(1,jj)=YpBar2(jj);
            YpontSplitBarDx(2,jj)=YpBar2(jj);
            

            XpontSplitSx(1,ii)=NaN;
            XpontSplitSx(2,ii)=NaN;
            YpontSplitSx(1,ii)=NaN;
            YpontSplitSx(2,ii)=NaN;
            XpontSplitDx(1,ii)=NaN;
            XpontSplitDx(2,ii)=NaN;
            YpontSplitDx(1,ii)=NaN;
            YpontSplitDx(2,ii)=NaN;
            
            
            temp.XpontRadDx(ii)    = XpontRadDx(ii);
            temp.YpontRadDx(ii)    = YpontRadDx(ii);
            temp.XpontRadSx(ii)    = XpontRadSx(ii);
            temp.YpontRadSx(ii)    = YpontRadSx(ii);
            temp.XpontRadBarDx(ii) = XpontRadBarDx(ii);
            temp.XpontRadBarSx(ii) = XpontRadBarSx(ii);
            temp.YpontRadBarDx(ii) = YpontRadBarDx(ii);
            temp.YpontRadBarSx(ii) = YpontRadBarSx(ii);
            temp.splitRib(ii)      = splitRib(ii);
            
            temp.XpontSplitBarSx(:,ii) = XpontSplitBarSx(:,ii);
            temp.YpontSplitBarSx(:,ii) = YpontSplitBarSx(:,ii);
            temp.XpontSplitBarDx(:,ii) = XpontSplitBarDx(:,ii);
            temp.YpontSplitBarDx(:,ii) = YpontSplitBarDx(:,ii);
            temp.XpontSplitDx(:,ii)    = XpontSplitDx(:,ii);
            temp.YpontSplitDx(:,ii)    = YpontSplitDx(:,ii);
            temp.XpontSplitSx(:,ii)    = XpontSplitSx(:,ii);
            temp.YpontSplitSx(:,ii)    = YpontSplitSx(:,ii);
            
            temp.xS01k(ii) = xS01k(ii);
            temp.yS01k(ii) = yS01k(ii);
            temp.xS02k(ii) = xS02k(ii);
            temp.yS02k(ii) = yS02k(ii);
            %
            temp.xI01k(ii) = xI01k(ii);
            temp.yI01k(ii) = yI01k(ii);
            temp.xI02k(ii) = xI02k(ii);
            temp.yI02k(ii) = yI02k(ii);
            
            temp.flag_cent(ii)    = flag_cent(ii);
            temp.flag_ext(ii)     = flag_ext(ii);
            temp.flag_shift(ii)   = flag_shift(ii);
            temp.flag_shiftUP(ii) = flag_shiftUP(ii);
            
            geo.pontRang(ii)=pontRang(ii);
            geo.pontRoffset(ii)=pontRoffset(ii);
            geo.RotorFillet(2,(ii))=RotorFilletRad1(ii);
            geo.RotorFillet(3,(ii))=RotorFilletRad2(ii);
            geo.pontR(ii)=pont(ii);
            geo.RotorFilletRad1(ii) = RotorFilletRad1(ii);
            geo.RotorFilletRad2(ii) = RotorFilletRad2(ii);
        end
        
        geo.pontR = pont;
        
        
        if any(radial_ribs_split==1)
            for jj=1:nlay
                if pont(jj)<pont0
                    pont(jj)=0;
                elseif pont(jj)<2*pont0
                    pont(jj)=2*pont0;
                end
            end
            
            for ii=1:nlay
                if pont(ii)>0
                    splitRib(ii)=1;
                end
                XpontRadBarSx(ii)=B1k(ii);
                YpontRadBarSx(ii)=0;
                XpontRadBarDx(ii)=B2k(ii);
                YpontRadBarDx(ii)=0;
                XpontRadDx(ii)=NaN;
                YpontRadDx(ii)=0;
                XpontRadSx(ii)=NaN;
                YpontRadSx(ii)=0;
            end
            
            %no double fillet radius in split ribs
            ii=(radial_ribs_split==1);
            geo.RotorFillet2(ii) = RotorFilletRad1(ii);
            
            index_ang=pontRang>(90/geo.p);
            if any(index_ang)
                disp('Limited Radial Ribs Angle @ 90/p')
                pontRang(index_ang)=(90/geo.p);
            end
            pontRang(pontRang<-85)=-85;
            m1=tand(pontRang);
            
            
            if any(pontRoffset)
                pontRoffset(pontRoffset<RotorFilletRad & (pontRoffset~=0))=RotorFilletRad(pontRoffset<RotorFilletRad & (pontRoffset~=0));
            end
            
            index_limOff=pontRoffset>0;
            if any(index_limOff)
                [a,b,c]=retta_per_2pti(XpBar1,YpBar1,XpBar2,YpBar2);
                [m_lim,q_lim,~]=retta_abc2mq(a,b,c);
                m1(m1>m_lim & index_limOff)=m_lim(m1>m_lim & index_limOff);
                pontRang=round(atand(m1),2);
                m1=tand(pontRang);
            end
            
            flag_cent=zeros(1,nlay);
            flag_ext=zeros(1,nlay);
            flag_shift=zeros(1,nlay);
            flag_shiftUP=zeros(1,nlay);
            
            
            %%%%cambiare pontRoff
            %pontRoffset=RotorFilletRad+pontRoffset;
            pontRoffset=-pontRoffset;
            %                 index_shift= StartRibs>0;
            %                 if sum(index_shift)
            
            %                 end
            
            index_shiftUP= pontRoffset<0;
            if any(index_shiftUP)
                flag_shiftUP(index_shiftUP)=1;
                %m1=tand(ang_rib);
                %q1=(YpBar2-StartRibs)-m1.*XpBar2;
                index_limang= pontRang>=0 & index_shiftUP;
                if any(index_limang)
                    pontRang(index_limang)=-0.1;
                    disp('Limited  angle at the Central Ribs (Lateral Ribs are feasible only with negative angle)')
                end
                
            end
            dist1=((xxD2k - XpBar2).^2 + (yyD2k - YpBar2).^2).^0.5-RotorFilletRad-pont/2;
            index_lim= abs(pontRoffset)>dist1;
            if any(index_lim)
                pontRoffset(index_lim & pontRoffset<0)=-round(0.9*abs(dist1(index_lim & pontRoffset<0)),2);
                pontRoffset(index_lim & pontRoffset>0)=round(0.9*abs(dist1(index_lim & pontRoffset>0)),2);
                %disp('Limited radial ribs offset (#1)')
                index_lim0=dist1(index_lim)<0;
                if any(index_lim0)
                    pontRoffset(index_lim0)=0;
                    flag_shiftUP(index_lim0)=0;
                end
            end
            flag_shift(pontRoffset>0)=1;
            
            
            q1=(YpBar2-pontRoffset)-m1.*XpBar2;
            
            
            
            index=m1.*XpBar1+q1<YpBar1 &not(flag_shiftUP);
            if any(index)
                flag_cent(index)=1;
                flag_ext(index)=0;
                
                index_up = flag_shiftUP==1 & index==1;
                %                     index_up1= YpontSplitSx(2,(index_up))<YpBar1(index_up);
                %                         if sum(index_up1)
                %                             disp('The selected radial ribs offset creates interference and it should be increased')
                %                             pontRoffset(index_up1)=0;
                %                             %                             xp2(index_up1)=XpBar1(index_up1);
                %                             %                             yp2(index_up1)=YpBar1(index_up1);
                %                             %
                %                             %                             q2dx=YpBar2-m2.*XpBar2;
                %                             %                             %q24(index_up)=yp4(index_up)-m1(index_up).*xp4(index_up);
                %                             %                             xp4(index_up1)=(q2dx(index_up1)-q24(index_up1))./(m1(index_up1)-m2(index_up1));
                %                             %                             yp4(index_up1)=m1(index_up1).*xp4(index_up1)+q24(index_up1);
                %                             %
                %                             %                              a(index_up1)=(1+m2(index_up1).*m2(index_up1));
                %                             %                              b(index_up1)=-2*xp4(index_up1)+2*m2(index_up1).*(q4(index_up1)-yp4(index_up1));
                %                             %                              c(index_up1)=-pont(index_up1)/2.*pont(index_up1)/2+xp4(index_up1).^2+(q4(index_up1)-yp4(index_up1)).^2;
                %                             %                              xp3(index_up1)=(-b(index_up1)+sqrt(b(index_up1).^2-4*a(index_up1).*c(index_up1)))./(2*a(index_up1));
                %                             %                              yp3(index_up1)=m2(index_up1).*xp3(index_up1)+q4(index_up1);
                %                         end
                
                
                xp3(index)=XpBar2(index);
                yp3(index)=YpBar2(index)-pontRoffset(index);
                
                xp1(index)=XpBar1(index);
                yp1(index)=m1(index).*xp1(index)+q1(index);
                
                %                     index_check=yp1(index)>YpBar1(index);
                %                     if any(index_check)~=0
                %                     [a(index_check),b(index_check),c(index_check)]=retta_per_2pti(XpBar1(index_check),YpBar1(index_check),xxD1k(index_check),yyD1k(index_check));
                %                     [m2(index_check),q2(index_check),~]=retta_abc2mq(a(index_check),b(index_check),c(index_check));
                %
                %                     xp1(index_check)=(q2(index_check)-q1(index_check))./(m1(index_check)-m2(index_check));
                %                     yp1(index_check)=m2(index_check).*xp1(index_check)+q2(index_check);
                %
                %                     flag_ext(index_check)=1;
                %
                %                     end
                
                
                
                dist((index))=((xp1(index) - XpBar1(index)).^2 + (yp1(index) - YpBar1(index)).^2).^0.5;
                index1= (RotorFilletRad(index)>dist(index));
                if any(index1)
                    flag_ext(index1)=1;
                    index(index1)=0;
                    %                         disp('Limited radial ribs fillet (#0)')
                    %                         RotorFilletRad(index1)=dist(index1);
                end
                
                xp2(index)=xp1(index);
                yp2(index)=yp1(index)-pont(index)/2;
                
                
                
                xp4(index)=xp3(index);
                yp4(index)=yp3(index)-pont(index)/2;
                
                %fillet
                %rotor fillet alla barriera laterale
                %                     [a,b,c]=retta_per_2pti(XpBar2,YpBar2,xxD2k,yyD2k);
                %                     [m2,q2,~]=retta_abc2mq(a,b,c);  %buttare ~
                %                     [a(index),b(index),c(index)]=retta_per_2pti(XpBar1(index),YpBar1(index),xxD1k(index),yyD1k(index));
                %                     [m2(index),q2(index),~]=retta_abc2mq(a(index),b(index),c(index));  %buttare ~
                %                     q2(index)=YpBar2(index)-m2(index).*XpBar2(index);
                
                [a,b,c]=retta_per_2pti(XpBar1,YpBar1,xxD1k,yyD1k);
                [m2,q2,~]=retta_abc2mq(a,b,c);
                q2=YpBar2-m2.*XpBar2;
                
                %                     array1=sum(index);
                %                     m1perp(index)=-ones(1,nlay)./m1(index);
                %                     m2perp(index)=-ones(1,nlay)./m2(index);
                m1perp=-ones(1,nlay)./m1;
                m2perp=-ones(1,nlay)./m2;
                
                
                
                %                     a=(1+m2.*m2);
                %                     b=-2*xp1+2*m2.*(q2-yp1);
                %                     c=-RotorFilletRad.*RotorFilletRad+xp1.^2+(q2-yp1).^2;
                XpontSplitBarSx(1,(index))=xp1(index);
                YpontSplitBarSx(1,(index))=yp1(index)+RotorFilletRad(index);
                
                %                     if m1(index)~=0
                a(index)=(1+m1(index).*m1(index));
                b(index)=-2*xp1(index)+2*m1(index).*(q1(index)-yp1(index));
                c(index)=-RotorFilletRad(index).*RotorFilletRad(index)+xp1(index).^2+(q1(index)-yp1(index)).^2;
                XpontSplitSx(1,(index))=(-b(index)+sqrt(b(index).^2-4*a(index).*c(index)))./(2*a(index));
                YpontSplitSx(1,(index))=m1(index).*XpontSplitSx(1,(index))+q1(index);
                %                     else
                %                         XpontSplitSx(1,(index))=xp1(index)+RotorFilletRad(index);
                %                         YpontSplitSx(1,(index))=yp1(index);
                %                     end
                
                index_st= flag_shift & index;
                if any(index_st)
                    %                         a(index_st)=(1+m2(index_st)).*m2(index_st);
                    %                         b(index_st)=-2*xp3(index_st)+2*m2(index_st).*(q2(index_st))-yp3(index_st);
                    %                         c(index_st)=-RotorFilletRad(index_st).*RotorFilletRad(index_st)+xp3(index_st).^2+(q2(index)-yp3(index_st)).^2;
                    %                         XpontSplitBarDx(1,(index_st))=(-b(index_st))+sqrt(b.^2-4*a(index_st)).*c(index_st)./(2*a(index_st));
                    %                         YpontSplitBarDx(1,(index_st))=m2(index_st).*XpontSplitBarDx(1,(index_st))+q2(index_st);
                    XpontSplitBarDx(1,(index_st))=xp3(index_st);
                    YpontSplitBarDx(1,(index_st))=yp3(index_st)+RotorFilletRad(index_st);
                    %
                    %                         a(index_st)=(1+m1(index_st).*m1(index_st);
                    %                         b(index_st)=-2*xp3(index_st)+2*m1(index_st).*(q1(index_st)-yp3(index_st));
                    %                         c(index_st)=-RotorFilletRad(index_st).*RotorFilletRad(index_st)+xp3(index_st).^2+(q1(index_st)-yp3(index_st)).^2;
                    %                         XpontSplitDx(1,(index_st))=(-b(index_st)-sqrt(b(index_st).^2-4*a(index_st)).*c(index_st))./(2*a(index_st));
                    %                         YpontSplitDx(1,(index_st))=m1(index_st).*XpontSplitDx(1,(index_st))+q1(index_st);
                    %                         if m1(index_st)~=0
                    a(index_st)=(1+m1(index_st).*m1(index_st));
                    b(index_st)=-2*xp3(index_st)+2*m1(index_st).*(q1(index_st)-yp3(index_st));
                    c(index_st)=-RotorFilletRad(index_st).*RotorFilletRad(index_st)+xp3(index_st).^2+(q1(index_st)-yp3(index_st)).^2;
                    XpontSplitDx(1,(index_st))=(-b(index_st)-sqrt(b(index_st).^2-4*a(index_st).*c(index_st)))./(2*a(index_st));
                    YpontSplitDx(1,(index_st))=m1(index_st).*XpontSplitDx(1,(index_st))+q1(index_st);
                    %                         else
                    %                             XpontSplitDx(1,(index_st))=xp3(index_st)-RotorFilletRad(index_st);
                    %                             YpontSplitDx(1,(index_st))=yp3(index_st);
                    %                         end
                end
                
                flag_shift_n=not(flag_shift);
                index_st1= flag_shift_n & index;
                if any(index_st1)
                    %                         XpontSplitBarDx(1,(index_st1))=xp3(index_st1);
                    %                         YpontSplitBarDx(1,(index_st1))=yp3(index_st1)+RotorFilletRad(index_st1);
                    
                    a(index_st1)=(1+m2(index_st1).*m2(index_st1));
                    b(index_st1)=-2*xp3(index_st1)+2*m2(index_st1).*(q2(index_st1)-yp3(index_st1));
                    c(index_st1)=-RotorFilletRad(index_st1).*RotorFilletRad(index_st1)+xp3(index_st1).^2+(q2(index_st1)-yp3(index_st1)).^2;
                    XpontSplitBarDx(1,(index_st1))=(-b(index_st1)+sqrt(b(index_st1).^2-4*a(index_st1).*c(index_st1)))./(2*a(index_st1));
                    YpontSplitBarDx(1,(index_st1))=m2(index_st1).*XpontSplitBarDx(1,(index_st1))+q2(index_st1);
                    
                    
                    %                         if m1(index_st1)~=0
                    a(index_st1)=(1+m1(index_st1).*m1(index_st1));
                    b(index_st1)=-2*xp3(index_st1)+2*m1(index_st1).*(q1(index_st1)-yp3(index_st1));
                    c(index_st1)=-RotorFilletRad(index_st1).*RotorFilletRad(index_st1)+xp3(index_st1).^2+(q1(index_st1)-yp3(index_st1)).^2;
                    XpontSplitDx(1,(index_st1))=(-b(index_st1)-sqrt(b(index_st1).^2-4*a(index_st1).*c(index_st1)))./(2*a(index_st1));
                    YpontSplitDx(1,(index_st1))=m1(index_st1).*XpontSplitDx(1,(index_st1))+q1(index_st1);
                    %                         else
                    %                             XpontSplitDx(1,(index_st1))=xp3(index_st1)-RotorFilletRad(index_st1);
                    %                             YpontSplitDx(1,(index_st1))=yp3(index_st1);
                    %                         end
                end
                
                %index_inf=zeros(1,nlay);
                index_inf=(isfinite(m1perp));
                if sum(index_inf)>0
                    q22(index_inf)=YpontSplitSx(1,(index_inf))-m1perp(index_inf).*XpontSplitSx(1,(index_inf));
                    yS01k(index_inf)=YpontSplitBarSx(1,(index_inf));
                    xS01k(index_inf)=(yS01k(index_inf)-q22(index_inf))./m1perp(index_inf);
                    
                    index_inf_s = (pontRoffset==0 & index_inf);
                    if sum(index_inf_s)>0
                        q33=zeros(1,nlay);
                        q44=zeros(1,nlay);
                        q33(index_inf_s)=YpontSplitBarDx(1,(index_inf_s))-m2perp(index_inf_s).*XpontSplitBarDx(1,(index_inf_s));
                        q44(index_inf_s)=YpontSplitDx(1,(index_inf_s))-m1perp(index_inf_s).*XpontSplitDx(1,(index_inf_s));
                        xS02k(index_inf_s)=(q44(index_inf_s)-q33(index_inf_s))./(m2perp(index_inf_s)-m1perp(index_inf_s));
                        yS02k(index_inf_s)=m2perp(index_inf_s).*xS02k(index_inf_s)+q33(index_inf_s);
                    end
                    index_inf_sn = (not(index_inf_s) & index_inf);
                    if sum(index_inf_sn)>0
                        q44(index_inf_sn)=YpontSplitDx(1,(index_inf_sn))-m1perp(index_inf_sn).*XpontSplitDx(1,(index_inf_sn));
                        yS02k(index_inf_sn)=YpontSplitBarDx(1,(index_inf_sn));
                        xS02k(index_inf_sn)=(yS02k(index_inf_sn)-q44(index_inf_sn))./m1perp(index_inf_sn);
                    end
                end
                
                index_infn=not(index_inf);
                if sum(index_infn)>0
                    xS01k(index_infn)=XpontSplitSx(1,(index_infn));
                    yS01k(index_infn)=YpontSplitBarSx(1,(index_infn));
                    
                    index_inf1_s = (pontRoffset==0 & index_infn);
                    if sum(index_inf1_s)>0
                        q33(index_inf1_s)=YpontSplitBarDx(1,(index_inf1_s))-m2perp(index_inf1_s).*XpontSplitBarDx(1,(index_inf1_s));
                        xS02k(index_inf1_s)=XpontSplitDx(1,(index_inf1_s));
                        yS02k(index_inf1_s)=m2perp(index_inf1_s).*xS02k(index_inf1_s)+q33(index_inf1_s);
                    end
                    
                    index_inf1_sn=(not(index_inf1_s) & index_infn);
                    if sum(index_inf1_sn)>0
                        xS02k(index_inf1_sn)=XpontSplitDx(1,(index_inf1_sn));
                        yS02k(index_inf1_sn)=YpontSplitBarDx(1,(index_inf1_sn));
                    end
                end
                
                q11(index)=yp2(index)-m1(index).*xp2(index);
                
                a(index)=(1+m1(index).*m1(index));
                b(index)=-2*xp2(index)+2*m1(index).*(q11(index)-yp2(index));
                c(index)=-RotorFilletRad(index).*RotorFilletRad(index)+xp2(index).^2+(q11(index)-yp2(index)).^2;
                XpontSplitSx(2,(index))=(-b(index)+sqrt(b(index).^2-4*a(index).*c(index)))./(2*a(index));
                YpontSplitSx(2,(index))=m1(index).*XpontSplitSx(2,(index))+q11(index);
                
                a(index)=(1+m1(index).*m1(index));
                b(index)=-2*xp4(index)+2*m1(index).*(q11(index)-yp4(index));
                c(index)=-RotorFilletRad(index).*RotorFilletRad(index)+xp4(index).^2+(q11(index)-yp4(index)).^2;
                XpontSplitDx(2,(index))=(-b(index)-sqrt(b(index).^2-4*a(index).*c(index)))./(2*a(index));
                YpontSplitDx(2,(index))=m1(index).*XpontSplitDx(2,(index))+q11(index);
                
                XpontSplitBarSx(2,(index))=xp2(index);
                YpontSplitBarSx(2,(index))=yp2(index)-RotorFilletRad(index);
                
                XpontSplitBarDx(2,(index))=xp4(index);
                YpontSplitBarDx(2,(index))=yp4(index)-RotorFilletRad(index);
                
                %index3=zeros(1,nlay);
                index3= isfinite(m1perp) & index;
                
                if sum(index3)>0
                    q22(index3)=YpontSplitSx(2,(index3))-m1perp(index3).*XpontSplitSx(2,(index3));
                    yI01k(index3)=YpontSplitBarSx(2,(index3));
                    xI01k(index3)=(yI01k(index3)-q22(index3))./m1perp(index3);
                    
                    q44(index3)=YpontSplitDx(2,(index3))-m1perp(index3).*XpontSplitDx(2,(index3));
                    yI02k(index3)=YpontSplitBarDx(2,(index3));
                    xI02k(index3)=(yI02k(index3)-q44(index3))./m1perp(index3);
                    
                end
                
                index3=not(index3);
                if sum(index3)>0
                    xI01k(index3)=XpontSplitSx(2,(index3));
                    yI01k(index3)=YpontSplitBarSx(2,(index3));
                    
                    xI02k(index3)=XpontSplitDx(2,(index3));
                    yI02k(index3)=YpontSplitBarDx(2,(index3));
                end
                
                
                %index2=zeros(1,nlay);
                index2= (pontRang(index)>0 & flag_shift(index)==0);
                if sum(index2)>0
                    disp('Limited radial ribs fillet (#1)')
                    xS02k(index2)=NaN;
                    yS02k(index2)=NaN;
                    XpontSplitBarDx(1,index2)=XpBar2(index2);
                    YpontSplitBarDx(1,index2)=YpBar2(index2);
                    
                    XpontSplitDx(1,index2)=XpBar2(index2);
                    YpontSplitDx(1,index2)=YpBar2(index2);
                end
            end
            
            index=not(index);
            if any(index)
                
                flag_ext(index)=1;
                [a(index),b(index),c(index)]=retta_per_2pti(XpBar1(index),YpBar1(index),xxD1k(index),yyD1k(index));
                [m2(index),q2(index),~]=retta_abc2mq(a(index),b(index),c(index));
                
                
                xp3(index)=XpBar2(index);
                yp3(index)=YpBar2(index)-pontRoffset(index);
                
                xp1(index)=(q2(index)-q1(index))./(m1(index)-m2(index));
                yp1(index)=m1(index).*xp1(index)+q1(index);
                
                
                a(index)=(1+m2(index).*m2(index));
                b(index)=-2*xp1(index)+2*m2(index).*(q2(index)-yp1(index));
                c(index)=-pont(index)/2.*pont(index)/2+xp1(index).^2+(q2(index)-yp1(index)).^2;
                xp2(index)=(-b(index)-sqrt(b(index).^2-4*a(index).*c(index)))./(2*a(index));
                yp2(index)=m2(index).*xp2(index)+q2(index);
                q4(index)=yp2(index)-m1(index).*xp2(index);
                
                %                     xp2=(q2-q4)./(m1-m2);
                %                     yp2=m1.*xp2+q4;
                xp4(index)=xp3(index);
                yp4(index)=m1(index).*xp3(index)+q4(index);
                
                %                     xp4=XpBar2;
                %                     yp4=YpBar2-pont/2;
                
                
                index_up = flag_shiftUP==1 & index==1;
                if any(index_up)
                    q4(index_up)=YpBar2(index_up)-m2(index_up).*XpBar2(index_up);
                    
                    a(index_up)=(1+m2(index_up).*m2(index_up));
                    b(index_up)=-2*XpBar2(index_up)+2*m2(index_up).*(q4(index_up)-YpBar2(index_up));
                    c(index_up)=-pontRoffset(index_up).*pontRoffset(index_up)+XpBar2(index_up).^2+(q4(index_up)-YpBar2(index_up)).^2;
                    xp4(index_up)=(-b(index_up)+sqrt(b(index_up).^2-4*a(index_up).*c(index_up)))./(2*a(index_up));
                    yp4(index_up)=m2(index_up).*xp4(index_up)+q4(index_up);
                    
                    a(index_up)=(1+m2(index_up).*m2(index_up));
                    b(index_up)=-2*xp4(index_up)+2*m2(index_up).*(q4(index_up)-yp4(index_up));
                    c(index_up)=-pont(index_up)/2.*pont(index_up)/2+xp4(index_up).^2+(q4(index_up)-yp4(index_up)).^2;
                    xp3(index_up)=(-b(index_up)+sqrt(b(index_up).^2-4*a(index_up).*c(index_up)))./(2*a(index_up));
                    yp3(index_up)=m2(index_up).*xp3(index_up)+q4(index_up);
                    
                    q24(index_up)=yp4(index_up)-m1(index_up).*xp4(index_up);
                    xp2(index_up)=(q2(index_up)-q24(index_up))./(m1(index_up)-m2(index_up));
                    yp2(index_up)=m1(index_up).*xp2(index_up)+q24(index_up);
                    
                    %                         a(index_up)=(1+m1(index_up).*m1(index_up));
                    %                         b(index_up)=-2*xp2(index_up)+2*m1(index_up).*(q11(index_up)-yp2(index_up));
                    %                         c(index_up)=-RotorFilletRad(index_up).*RotorFilletRad(index_up)+xp2(index_up).^2+(q11(index_up)-yp2(index_up)).^2;
                    %                         XpontSplitSx(2,(index_up))=(-b(index_up)+sqrt(b(index_up).^2-4*a(index_up).*c(index_up)))./(2*a(index_up));
                    %                         YpontSplitSx(2,(index_up))=m1(index_up).*XpontSplitSx(2,(index_up))+q11(index_up);
                    
                    
                    q31(index_up)=yp3(index_up)-m1(index_up).*xp3(index_up);
                    q1(index_up)=q31(index_up);
                    xp1(index_up)=(q2(index_up)-q31(index_up))./(m1(index_up)-m2(index_up));
                    yp1(index_up)=m1(index_up).*xp1(index_up)+q31(index_up);
                    
                    
                end
                
                
                
                
                
                index1= (xp2<XpBar1) & not(flag_shiftUP);
                if any(index1)
                    %xp2(index1)<XpBar1(index1);
                    xp2(index1)=XpBar1(index1);
                    yp2(index1)=m1(index1).*xp2(index1)+q4(index1);
                    flag_cent(index1)=1;
                else
                    dist(index1)=((xp2(index1) - XpBar1(index1)).^2 + (yp2(index1) - YpBar1(index1)).^2).^0.5;
                    index1=(RotorFilletRad(1,(index1))>dist(index1));
                    if index1
                        disp('Limited radial ribs fillet (#2)')
                        RotorFilletRad(index1)=round(dist(index1)*0.9,2);
                    end
                end
                
                a(index)=(1+m2(index).*m2(index));
                b(index)=-2*xp1(index)+2*m2(index).*(q2(index)-yp1(index));
                c(index)=-RotorFilletRad(index).*RotorFilletRad(index)+xp1(index).^2+(q2(index)-yp1(index)).^2;
                XpontSplitBarSx(1,(index))=(-b(index)+sqrt(b(index).^2-4*a(index).*c(index)))./(2*a(index));
                YpontSplitBarSx(1,(index))=m2(index).*XpontSplitBarSx(1,(index))+q2(index);
                
                a(index)=(1+m1(index).*m1(index));
                b(index)=-2*xp1(index)+2*m1(index).*(q1(index)-yp1(index));
                c(index)=-RotorFilletRad(index).*RotorFilletRad(index)+xp1(index).^2+(q1(index)-yp1(index)).^2;
                XpontSplitSx(1,(index))=(-b(index)+sqrt(b(index).^2-4*a(index).*c(index)))./(2*a(index));
                YpontSplitSx(1,(index))=m1(index).*XpontSplitSx(1,(index))+q1(index);
                
                index_check=index & flag_shift & not(flag_shiftUP);
                if any(index_check)
                    XpontSplitBarDx(1,(index_check))=xp3(index_check);
                    YpontSplitBarDx(1,(index_check))=yp3(index_check)+RotorFilletRad(index_check);
                end
                
                
                flag_shift_n=not(flag_shift);
                index_check=index & flag_shift_n;
                if any(index_check)
                    q3(index_check)=-m2(index_check).*XpBar2(index_check)+YpBar2(index_check);
                    a(index_check)=(1+m2(index_check).*m2(index_check));
                    b(index_check)=-2*xp3(index_check)+2*m2(index_check).*(q3(index_check)-yp3(index_check));
                    c(index_check)=-RotorFilletRad(index_check).*RotorFilletRad(index_check)+xp3(index_check).^2+(q3(index_check)-yp3(index_check)).^2;
                    XpontSplitBarDx(1,(index_check))=(-b(index_check)+sqrt(b(index_check).^2-4*a(index_check).*c(index_check)))./(2*a(index_check));
                    YpontSplitBarDx(1,(index_check))=m2(index_check).*XpontSplitBarDx(1,(index_check))+q3(index_check);
                end
                
                a(index)=(1+m1(index).*m1(index));
                b(index)=-2*xp3(index)+2*m1(index).*(q1(index)-yp3(index));
                c(index)=-RotorFilletRad(index).*RotorFilletRad(index)+xp3(index).^2+(q1(index)-yp3(index)).^2;
                XpontSplitDx(1,(index))=(-b(index)-sqrt(b(index).^2-4*a(index).*c(index)))./(2*a(index));
                YpontSplitDx(1,(index))=m1(index).*XpontSplitDx(1,(index))+q1(index);
                
                array1=sum(index);
                m1perp(index)=-ones(1,array1)./m1(index);
                m2perp(index)=-ones(1,array1)./m2(index);
                
                index_check=index & flag_shift;
                if any(index_check)
                    yS02k(index_check)=YpontSplitBarDx(1,(index_check));
                    q44(index_check)=YpontSplitDx(1,(index_check))-m1perp(index_check).*XpontSplitDx(1,(index_check));
                    xS02k(index_check)=(yS02k(index_check)-q44(index_check))./m1perp(index_check);
                end
                
                index_check=index & flag_shift_n;
                if any(index_check)
                    q33(index_check)=YpontSplitBarDx(1,(index_check))-m2perp(index_check).*XpontSplitBarDx(1,(index_check));
                    q44(index_check)=YpontSplitDx(1,(index_check))-m1perp(index_check).*XpontSplitDx(1,(index_check));
                    xS02k(index_check)=(q44(index_check)-q33(index_check))./(m2perp(index_check)-m1perp(index_check));
                    yS02k(index_check)=m2perp(index_check).*xS02k(index_check)+q33(index_check);
                end
                
                xS02k(~isfinite(m1perp) & index)=XpontSplitBarDx(1,~isfinite(m1perp) & index);
                q33=YpontSplitBarDx(1,:)-m2perp.*XpontSplitBarDx(1,:);
                yS02k(~isfinite(m1perp) & index)=m2perp(~isfinite(m1perp) & index).*xS02k(isfinite(m1perp)==0 & index)+q33(isfinite(m1perp)==0 & index);
                
                
                q11(index)=YpontSplitBarSx(1,(index))-m2perp(index).*XpontSplitBarSx(1,(index));
                q22(index)=YpontSplitSx(1,(index))-m1perp(index).*XpontSplitSx(1,(index));
                xS01k(index)=(q22(index)-q11(index))./(m2perp(index)-m1perp(index));
                xS01k(~isfinite(m1perp) & index)=XpontSplitBarSx(1,~isfinite(m1perp) & index);
                yS01k(index)=m2perp(index).*xS01k(index)+q11(index);
                
                
                q11(index)=yp2(index)-m1(index).*xp2(index);
                
                if flag_cent(index)==0
                    dist(index)=((xp2(index) - XpBar1(index)).^2 + (yp2(index) - YpBar1(index)).^2).^0.5;
                    index1=(RotorFilletRad(1,:)>dist & index);
                    if any(index1)
                        disp('Limited radial ribs fillet (#3)')
                        RotorFilletRad(index1)=round(dist(index1)*0.9,2);
                    end
                end
                
                a(index)=(1+m1(index).*m1(index));
                b(index)=-2*xp2(index)+2*m1(index).*(q11(index)-yp2(index));
                c(index)=-RotorFilletRad(index).*RotorFilletRad(index)+xp2(index).^2+(q11(index)-yp2(index)).^2;
                XpontSplitSx(2,(index))=(-b(index)+sqrt(b(index).^2-4*a(index).*c(index)))./(2*a(index));
                YpontSplitSx(2,(index))=m1(index).*XpontSplitSx(2,(index))+q11(index);
                
                %index_up1= YpontSplitSx(2,(index_up))<YpBar1(index_up) & flag_shiftUP==1;
                if sum(YpontSplitSx(2,(index))<YpBar1(index) & flag_shiftUP(index)==1)
                    disp('The selected radial ribs offset creates interference and it should be increased #1')
                    %                             xp2(index_up1)=XpBar1(index_up1);
                    %                             yp2(index_up1)=YpBar1(index_up1);
                    %
                    %                             q2dx=YpBar2-m2.*XpBar2;
                    %                             %q24(index_up)=yp4(index_up)-m1(index_up).*xp4(index_up);
                    %                             xp4(index_up1)=(q2dx(index_up1)-q24(index_up1))./(m1(index_up1)-m2(index_up1));
                    %                             yp4(index_up1)=m1(index_up1).*xp4(index_up1)+q24(index_up1);
                    %
                    %                              a(index_up1)=(1+m2(index_up1).*m2(index_up1));
                    %                              b(index_up1)=-2*xp4(index_up1)+2*m2(index_up1).*(q4(index_up1)-yp4(index_up1));
                    %                              c(index_up1)=-pont(index_up1)/2.*pont(index_up1)/2+xp4(index_up1).^2+(q4(index_up1)-yp4(index_up1)).^2;
                    %                              xp3(index_up1)=(-b(index_up1)+sqrt(b(index_up1).^2-4*a(index_up1).*c(index_up1)))./(2*a(index_up1));
                    %                              yp3(index_up1)=m2(index_up1).*xp3(index_up1)+q4(index_up1);
                end
                
                
                a(index)=(1+m1(index).*m1(index));
                b(index)=-2*xp4(index)+2*m1(index).*(q11(index)-yp4(index));
                c(index)=-RotorFilletRad(index).*RotorFilletRad(index)+xp4(index).^2+(q11(index)-yp4(index)).^2;
                XpontSplitDx(2,(index))=(-b(index)-sqrt(b(index).^2-4*a(index).*c(index)))./(2*a(index));
                YpontSplitDx(2,(index))=m1(index).*XpontSplitDx(2,(index))+q11(index);
                
                XpontSplitBarDx(2,(index))=xp4(index);
                YpontSplitBarDx(2,(index))=yp4(index)-RotorFilletRad(index);
                
                %index_up=zeros(1,nlay);
                index_up= flag_shiftUP==1;
                if any(index_up)
                    q3(index_up)=-m2(index_up).*XpBar2(index_up)+YpBar2(index_up);
                    a(index_up)=(1+m2(index_up).*m2(index_up));
                    b(index_up)=-2*xp4(index_up)+2*m2(index_up).*(q3(index_up)-yp4(index_up));
                    c(index_up)=-RotorFilletRad(index_up).*RotorFilletRad(index_up)+xp4(index_up).^2+(q3(index_up)-yp4(index_up)).^2;
                    XpontSplitBarDx(2,(index_up))=(-b(index_up)-sqrt(b(index_up).^2-4*a(index_up).*c(index_up)))./(2*a(index_up));
                    YpontSplitBarDx(2,(index_up))=m2(index_up).*XpontSplitBarDx(2,(index_up))+q3(index_up);
                end
                
                index_c=flag_cent & flag_ext;
                if any(index_c)
                    XpontSplitBarSx(2,(index_c))=xp2(index_c);
                    YpontSplitBarSx(2,(index_c))=yp2(index_c)-RotorFilletRad(index_c);
                    
                    yI01k(index_c)=YpontSplitBarSx(2,(index_c));
                    q44(index_c)=YpontSplitSx(2,(index_c))-m1perp(index_c).*XpontSplitSx(2,(index_c));
                    xI01k(index_c)=(yI01k(index_c)-q44(index_c))./m1perp(index_c);
                end
                index_c=not(index_c) & flag_ext;
                if any(index_c)
                    a(index_c)=(1+m2(index_c).*m2(index_c));
                    b(index_c)=-2*xp2(index_c)+2*m2(index_c).*(q2(index_c)-yp2(index_c));
                    c(index_c)=-RotorFilletRad(index_c).*RotorFilletRad(index_c)+xp2(index_c).^2+(q2(index_c)-yp2(index_c)).^2;
                    XpontSplitBarSx(2,(index_c))=(-b(index_c)-sqrt(b(index_c).^2-4*a(index_c).*c(index_c)))./(2*a(index_c));
                    YpontSplitBarSx(2,(index_c))=m2(index_c).*XpontSplitBarSx(2,(index_c))+q2(index_c);
                    
                    q33(index_c)=YpontSplitBarSx(2,(index_c))-m2perp(index_c).*XpontSplitBarSx(2,(index_c));
                    q44(index_c)=YpontSplitSx(2,(index_c))-m1perp(index_c).*XpontSplitSx(2,(index_c));
                    xI01k(index_c)=(q44(index_c)-q33(index_c))./(m2perp(index_c)-m1perp(index_c));
                    yI01k(index_c)=m2perp(index_c).*xI01k(index_c)+q33(index_c);
                end
                
                index_up_i= index_up & index;
                if any(index_up_i)
                    q33(index_up_i)=YpontSplitBarDx(2,(index_up_i))-m2perp(index_up_i).*XpontSplitBarDx(2,(index_up_i));
                    q44(index_up_i)=YpontSplitDx(2,(index_up_i))-m1perp(index_up_i).*XpontSplitDx(2,(index_up_i));
                    xI02k(index_up_i)=(q44(index_up_i)-q33(index_up_i))./(m2perp(index_up_i)-m1perp(index_up_i));
                    yI02k(index_up_i)=m2perp(index_up_i).*xI02k(index_up_i)+q33(index_up_i);
                end
                
                index_down= not(index_up) & index;
                if any(index_down)
                    yI02k(index_down)=YpontSplitBarDx(2,(index_down));
                    q44(index_down)=YpontSplitDx(2,(index_down))-m1perp(index_down).*XpontSplitDx(2,(index_down));
                    xI02k(index_down)=(yI02k(index_down)-q44(index_down))./m1perp(index_down);
                end
                
                
                index2= (pontRang>0 & index);
                if any(index2)
                    disp('Limited radial ribs fillet (#4)')
                    xS02k(index2)=NaN;
                    yS02k(index2)=NaN;
                    
                    XpontSplitBarDx(1,index2)=XpBar2(index2);
                    YpontSplitBarDx(1,index2)=YpBar2(index2);
                    
                    XpontSplitDx(1,index2)=XpBar2(index2);
                    YpontSplitDx(1,index2)=YpBar2(index2);
                end
                
                
            end
            
            %index=zeros(1,nlay);
            %             index= (temp.xxD2k<XpontSplitBarDx(1,:)) &(temp.yyD2k<YpontSplitBarDx(1,:));
            %             if sum(index)>0
            %                  error('DECREASE THE RADIAL ROTOR FILLET!')
            %                 disp('Inferior arc limited in rotor fillet');
            % %                 XpontSplitBarDx(index)=temp.xxD2k(index);
            % %                 YpontSplitBarDx(index)=temp.yyD2k(index);
            % %                 RotorFilletLimited=sqrt((XpontSplitBarDx(index)-xp3(index)).^2+(YpontSplitBarDx(index)-yp3(index)).^2);
            % %                 q1=yp3-m1.*xp3;
            % %                 a(index)=(1+m1(index).*m1(index));
            % %                 b(index)=-2*xp3(index)+2*m1(index).*(q1(index)-yp3(index));
            % %                 c(index)=-RotorFilletLimited.*RotorFilletLimited+xp3(index).^2+(q1(index)-yp3(index)).^2;
            % %                 XpontSplitDx(1,index)=(-b(index)-sqrt(b(index).^2-4*a(index).*c(index)))./(2*a(index));
            % %                 YpontSplitDx(1,index)=m1(index).*XpontSplitDx(index)+q1(index);
            %             end
            
            index = (XpontSplitBarSx(1,:)<XpontRadBarSx);
            if sum(index)
                error('DECREASE THE RADIAL ROTOR FILLET!')
            end
            
            %
            for ii=1:nlay
                if pont(ii)==0
                    YpontSplitBarSx(2,ii)=YpBar1(ii);
                    YpontSplitBarDx(2,ii)=YpBar2(ii);
                    YpontRadBarSx(ii)=0;
                    XpontRadBarDx(ii)=B2k(ii);
                    YpontRadBarDx(ii)=0;
                    XpontRadDx(ii)=NaN;
                    YpontRadDx(ii)=0;
                    XpontRadSx(ii)=NaN;
                    YpontRadSx(ii)=0;
                else
                    
                end
            end
            
            XpontSplitSx(1,pont==0)=NaN;
            XpontSplitSx(2,pont==0)=NaN;
            YpontSplitSx(1,pont==0)=NaN;
            YpontSplitSx(2,pont==0)=NaN;
            XpontSplitDx(1,pont==0)=NaN;
            XpontSplitDx(2,pont==0)=NaN;
            YpontSplitDx(1,pont==0)=NaN;
            YpontSplitDx(2,pont==0)=NaN;
            
            pontRoffset=-pontRoffset;
            
            
            %assign values pontsplit
            ii=(radial_ribs_split==1);
            temp.XpontRadDx(ii)    = XpontRadDx(ii);
            temp.YpontRadDx(ii)    = YpontRadDx(ii);
            temp.XpontRadSx(ii)    = XpontRadSx(ii);
            temp.YpontRadSx(ii)    = YpontRadSx(ii);
            temp.XpontRadBarDx(ii) = XpontRadBarDx(ii);
            temp.XpontRadBarSx(ii) = XpontRadBarSx(ii);
            temp.YpontRadBarDx(ii) = YpontRadBarDx(ii);
            temp.YpontRadBarSx(ii) = YpontRadBarSx(ii);
            temp.splitRib(ii)      = splitRib(ii);
            
            temp.XpontSplitBarSx(:,ii) = XpontSplitBarSx(:,ii);
            temp.YpontSplitBarSx(:,ii) = YpontSplitBarSx(:,ii);
            temp.XpontSplitBarDx(:,ii) = XpontSplitBarDx(:,ii);
            temp.YpontSplitBarDx(:,ii) = YpontSplitBarDx(:,ii);
            temp.XpontSplitDx(:,ii)    = XpontSplitDx(:,ii);
            temp.YpontSplitDx(:,ii)    = YpontSplitDx(:,ii);
            temp.XpontSplitSx(:,ii)    = XpontSplitSx(:,ii);
            temp.YpontSplitSx(:,ii)    = YpontSplitSx(:,ii);
            
            temp.xS01k(ii) = xS01k(ii);
            temp.yS01k(ii) = yS01k(ii);
            temp.xS02k(ii) = xS02k(ii);
            temp.yS02k(ii) = yS02k(ii);
            
            temp.xI01k(ii) = xI01k(ii);
            temp.yI01k(ii) = yI01k(ii);
            temp.xI02k(ii) = xI02k(ii);
            temp.yI02k(ii) = yI02k(ii);
            
            temp.flag_cent(ii)    = flag_cent(ii);
            temp.flag_ext(ii)     = flag_ext(ii);
            temp.flag_shift(ii)   = flag_shift(ii);
            temp.flag_shiftUP(ii) = flag_shiftUP(ii);
            
            geo.pontRang(ii)=pontRang(ii);
            geo.pontRoffset(ii)=pontRoffset(ii);
            geo.RotorFillet(2,(ii))=RotorFilletRad(ii);
            geo.RotorFillet(3,(ii))=NaN;
            geo.pontR(ii)=pont(ii);
            geo.RotorFillet1(ii) = RotorFilletRad1(ii);
            geo.RotorFillet2(ii) = 0;
        end
        
        
end

geo.radial_ribs_split = radial_ribs_split;

% end
% Points for radial ribs
% temp.XpontRadDx    = XpontRadDx;
% temp.YpontRadDx    = YpontRadDx;
% temp.XpontRadSx    = XpontRadSx;
% temp.YpontRadSx    = YpontRadSx;
% temp.XpontRadBarDx = XpontRadBarDx;
% temp.XpontRadBarSx = XpontRadBarSx;
% temp.YpontRadBarDx = YpontRadBarDx;
% temp.YpontRadBarSx = YpontRadBarSx;
% temp.splitRib      = splitRib;
%
% temp.XpontSplitBarSx = XpontSplitBarSx; % coordinates of split inner ribs
% temp.YpontSplitBarSx = YpontSplitBarSx;
% temp.XpontSplitBarDx = XpontSplitBarDx;
% temp.YpontSplitBarDx = YpontSplitBarDx;
% temp.XpontSplitDx    = XpontSplitDx;
% temp.YpontSplitDx    = YpontSplitDx;
% temp.XpontSplitSx    = XpontSplitSx;
% temp.YpontSplitSx    = YpontSplitSx;
%
% temp.xS01k = xS01k;
% temp.yS01k = yS01k;
% temp.xS02k = xS02k;
% temp.yS02k = yS02k;
%
% temp.xI01k = xI01k;
% temp.yI01k = yI01k;
% temp.xI02k = xI02k;
% temp.yI02k = yI02k;
%
% temp.flag_cent    = flag_cent;
% temp.flag_ext     = flag_ext;
% temp.flag_shift   = flag_shift;
% temp.flag_shiftUP = flag_shiftUP;
%
% geo.pontRang=pontRang;
% geo.pontRoffset=pontRoffset;
% geo.RotorFillet(2,:)=RotorFilletRad;
% geo.pontR=pont;

