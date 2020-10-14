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

B1k    = temp.B1k;
B2k    = temp.B2k;

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
        
        if geo.radial_ribs_eval == 0
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
        Afe=cumsum(A);
        
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
        
        if geo.radial_ribs_eval == 0
            pont = F_centrifuga/(sigma_max * l);    % mm
        else
            pont = geo.pontR;
        end
        
        %         XpontRadBarSx = zeros(1,nlay);
        %         YpontRadBarSx = zeros(1,nlay);
        %         XpontRadBarDx = zeros(1,nlay);
        %         YpontRadBarDx = zeros(1,nlay);
        %         XpontRadDx    = zeros(1,nlay);
        %         YpontRadDx    = zeros(1,nlay);
        %         XpontRadSx    = zeros(1,nlay);
        %         YpontRadSx    = zeros(1,nlay);
        
        XpontSplitBarSx = nan(2,nlay);
        YpontSplitBarSx = nan(2,nlay);
        XpontSplitBarDx = nan(2,nlay);
        YpontSplitBarDx = nan(2,nlay);
        XpontSplitDx    = nan(2,nlay);
        YpontSplitDx    = nan(2,nlay);
        XpontSplitSx    = nan(2,nlay);
        YpontSplitSx    = nan(2,nlay);
        
        
        if geo.radial_ribs_split
            % Controllo spessore ponticello:
            % - pont<pont0          --> no ponticello
            % - pont0<=pont<2*pont0 --> disegno 2 ponticelli spessi pont0
            % - pont >=2*pont0      --> disegno 2 ponticelli spessi pont/2
            for jj=1:nlay
                if pont(jj)<pont0
                    pont(jj)=0;
                elseif pont(jj)<2*pont0
                    pont(jj)=2*pont0;
                end
            end
            
            hpont=pont/2;
            
            % check feasability: outer side of the barrier shorter than the inner
            splitRib=ones(size(pont));
            for ii=1:nlay
                if (YpBar1(ii)<(YpBar2(ii)+racc_pont(ii)))
                    splitRib(ii)=0;
                end
            end
            
            for ii=1:nlay
                if (hpont(ii)>0)
                    if ~splitRib(ii) % error flag --> there is no space for the split rib, so the single rib is adopted
                        XpontSplitBarSx(1,ii)=B1k(ii);
                        XpontSplitBarSx(2,ii)=B1k(ii);
                        YpontSplitBarSx(1,ii)=YpBar1(ii);
                        YpontSplitBarSx(2,ii)=YpontSplitBarSx(1,ii);
                        XpontSplitBarDx(1,ii)=B2k(ii);
                        XpontSplitBarDx(2,ii)=B2k(ii);
                        YpontSplitBarDx(1,ii)=YpBar2(ii);
                        YpontSplitBarDx(2,ii)=YpontSplitBarDx(1,ii);
                        XpontSplitSx(1,ii)=NaN;
                        XpontSplitSx(2,ii)=NaN;
                        YpontSplitSx(1,ii)=NaN;
                        YpontSplitSx(2,ii)=NaN;
                        XpontSplitDx(1,ii)=NaN;
                        XpontSplitDx(2,ii)=NaN;
                        YpontSplitDx(1,ii)=NaN;
                        YpontSplitDx(2,ii)=NaN;
                        
                        XpontRadBarSx(ii)=B1k(ii);
                        YpontRadBarSx(ii)=hpont(ii)+racc_pont(ii);
                        XpontRadBarDx(ii)=B2k(ii);
                        YpontRadBarDx(ii)=hpont(ii)+racc_pont(ii);
                        XpontRadDx(ii)=B2k(ii)-racc_pont(ii);
                        YpontRadDx(ii)=hpont(ii);
                        XpontRadSx(ii)=B1k(ii)+racc_pont(ii);
                        YpontRadSx(ii)=hpont(ii);
                    else
                        XpontSplitBarSx(1,ii)=B1k(ii);
                        XpontSplitBarSx(2,ii)=B1k(ii);
                        YpontSplitBarSx(1,ii)=YpBar2(ii);
                        YpontSplitBarSx(2,ii)=YpontSplitBarSx(1,ii)-2*racc_pont(ii)-hpont(ii);
                        XpontSplitBarDx(1,ii)=B2k(ii);
                        XpontSplitBarDx(2,ii)=B2k(ii);
                        YpontSplitBarDx(1,ii)=YpBar2(ii);
                        YpontSplitBarDx(2,ii)=YpontSplitBarDx(1,ii)-2*racc_pont(ii)-hpont(ii);
                        XpontSplitSx(1,ii)=B1k(ii)+racc_pont(ii);
                        XpontSplitSx(2,ii)=B1k(ii)+racc_pont(ii);
                        YpontSplitSx(1,ii)=YpBar2(ii)-racc_pont(ii);
                        YpontSplitSx(2,ii)=YpontSplitSx(1,ii)-hpont(ii);
                        XpontSplitDx(1,ii)=B2k(ii)-racc_pont(ii);
                        XpontSplitDx(2,ii)=B2k(ii)-racc_pont(ii);
                        YpontSplitDx(1,ii)=YpBar2(ii)-racc_pont(ii);
                        YpontSplitDx(2,ii)=YpontSplitDx(1,ii)-hpont(ii);
                        
                        XpontRadBarSx(ii)=B1k(ii);
                        YpontRadBarSx(ii)=0;
                        XpontRadBarDx(ii)=B2k(ii);
                        YpontRadBarDx(ii)=0;
                        XpontRadDx(ii)=NaN;
                        YpontRadDx(ii)=0;
                        XpontRadSx(ii)=NaN;
                        YpontRadSx(ii)=0;
                        
                    end
                else
                    XpontSplitBarSx(1,ii)=B1k(ii);
                    XpontSplitBarSx(2,ii)=B1k(ii);
                    YpontSplitBarSx(1,ii)=YpBar2(ii);
                    YpontSplitBarSx(2,ii)=YpontSplitBarSx(1,ii);
                    XpontSplitBarDx(1,ii)=B2k(ii);
                    XpontSplitBarDx(2,ii)=B2k(ii);
                    YpontSplitBarDx(1,ii)=YpBar2(ii);
                    YpontSplitBarDx(2,ii)=YpontSplitBarDx(1,ii);
                    XpontSplitSx(1,ii)=NaN;
                    XpontSplitSx(2,ii)=NaN;
                    YpontSplitSx(1,ii)=NaN;
                    YpontSplitSx(2,ii)=NaN;
                    XpontSplitDx(1,ii)=NaN;
                    XpontSplitDx(2,ii)=NaN;
                    YpontSplitDx(1,ii)=NaN;
                    YpontSplitDx(2,ii)=NaN;
                    
                    XpontRadBarSx(ii)=B1k(ii);
                    YpontRadBarSx(ii)=0;
                    XpontRadBarDx(ii)=B2k(ii);
                    YpontRadBarDx(ii)=0;
                    XpontRadDx(ii)=NaN;
                    YpontRadDx(ii)=0;
                    XpontRadSx(ii)=NaN;
                    YpontRadSx(ii)=0;
                    
                end
            end
            
        else
            for jj=1:nlay
                if (pont(jj) < pont0) % non disegno i ponticelli radiali il cui spessore è minore della tolleranza di lavorazione per gli altri tipi di rotore
                    pont(jj)=0;
                end
            end
            
            hpont=pont/2;
            rac_pont=abs(B1k-B2k)/4;
            rac_pont=pont0;
            
            for ii=1:nlay
                if hpont(ii)>0
                    XpontRadBarSx(ii)=B1k(ii);
                    YpontRadBarSx(ii)=hpont(ii)+racc_pont(ii);
                    XpontRadBarDx(ii)=B2k(ii);
                    YpontRadBarDx(ii)=hpont(ii)+racc_pont(ii);
                    XpontRadDx(ii)=B2k(ii)-racc_pont(ii);
                    YpontRadDx(ii)=hpont(ii);
                    XpontRadSx(ii)=B1k(ii)+racc_pont(ii);
                    YpontRadSx(ii)=hpont(ii);
                else
                    XpontRadBarSx(ii)=B1k(ii);
                    YpontRadBarSx(ii)=0;
                    XpontRadBarDx(ii)=B2k(ii);
                    YpontRadBarDx(ii)=0;
                    XpontRadDx(ii)=NaN;
                    YpontRadDx(ii)=0;
                    XpontRadSx(ii)=NaN;
                    YpontRadSx(ii)=0;
                end
                
                XpontSplitBarSx(1,ii)=B1k(ii);
                XpontSplitBarSx(2,ii)=B1k(ii);
                YpontSplitBarSx(1,ii)=YpBar2(ii);
                YpontSplitBarSx(2,ii)=YpontSplitBarSx(1,ii);
                XpontSplitBarDx(1,ii)=B2k(ii);
                XpontSplitBarDx(2,ii)=B2k(ii);
                YpontSplitBarDx(1,ii)=YpBar2(ii);
                YpontSplitBarDx(2,ii)=YpontSplitBarDx(1,ii);
                XpontSplitSx(1,ii)=NaN;
                XpontSplitSx(2,ii)=NaN;
                YpontSplitSx(1,ii)=NaN;
                YpontSplitSx(2,ii)=NaN;
                XpontSplitDx(1,ii)=NaN;
                XpontSplitDx(2,ii)=NaN;
                YpontSplitDx(1,ii)=NaN;
                YpontSplitDx(2,ii)=NaN;
                
            end
            
        end
        temp.XpontSplitBarSx = XpontSplitBarSx; % coordinates of split inner ribs
        temp.YpontSplitBarSx = YpontSplitBarSx;
        temp.XpontSplitBarDx = XpontSplitBarDx;
        temp.YpontSplitBarDx = YpontSplitBarDx;
        temp.XpontSplitDx    = XpontSplitDx;
        temp.YpontSplitDx    = YpontSplitDx;
        temp.XpontSplitSx    = XpontSplitSx;
        temp.YpontSplitSx    = YpontSplitSx;
        
        geo.pontR = pont;
        
end

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


