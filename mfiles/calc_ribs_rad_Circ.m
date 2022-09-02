
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

function [temp,geo] = calc_ribs_rad_Circ(geo,mat,temp)

% NEW: function version of script calc_ribs_rad
% OK for CIRCULAR
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


RotorFilletRad = geo.RotorFillet1;
RotorFilletRad1 = geo.RotorFillet1;
RotorFilletRad2 = geo.RotorFillet2;
% delta_FBS       = geo.delta_FBS;

B1k    = temp.B1k;
B2k    = temp.B2k;

% xxD1k = temp.xxD1k;
% yyD1k = temp.yyD1k;
% xxD2k = temp.xxD2k;
yyD2k = temp.yyD2k;


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



RotorFilletRad2(RotorFilletRad2>hc/2)  = hc((RotorFilletRad2>hc/2))/2;
RotorFilletRad2(RotorFilletRad2<pont0) = pont0;
RotorFilletRad1(RotorFilletRad1>hc/2)  = hc((RotorFilletRad1>hc/2))/2;
RotorFilletRad1(RotorFilletRad1<pont0) = pont0;

%limit Radial ribs Fillet
RotorFilletRad(RotorFilletRad>hc/2)=round(hc((RotorFilletRad>hc/2))/2,2);
RotorFilletRad(RotorFilletRad<pont0)=pont0;

%check kOB with split ribs (angle feasible)


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
if ~isempty(yyD2k)
    pont(pont>(yyD2k-pont0)) = yyD2k(pont>(yyD2k-pont0))-pont0;
end

racc_pont = racc_pont*ones(1,nlay);
racc_pont(2*racc_pont>hc) = hc(2*racc_pont>hc)/3;
%         racc_pont = (RotorFilletRad2+RotorFilletRad1)/2;
%         RotorFilletRad2 = racc_pont;
%         RotorFilletRad1 = racc_pont;
racc_pont = mean([RotorFilletRad1;RotorFilletRad2]);
RotorFilletRad1 = racc_pont;
RotorFilletRad2 = racc_pont;

%%% Radial Ribs Node Coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jj = 1 : 2 : length(r_all)
    x10 = x0 - r_all(jj) - racc_pont(ceil(jj/2));   % (Ponticello radiale_vertice in alto a dx)=(Punto in basso a dx del raccordo))
    x11 = x0 - r_all(jj+1) + racc_pont(ceil(jj/2)); % (Ponticello radiale_vertice in alto a sx)=(Punto in basso a sx del raccordo))
    if pont(ceil(jj/2))>0
        if (x11 < x10)
            ii=ceil(jj/2);
            % Se i due punti non si incrociano (raccordo non troppo
            % grande rispetto alla larghezza barriera)
            hpont = pont(ii);

            [xS02k(ii),yS02k(ii)] = intersezione_retta_circonferenza(x0,0,-(r_all(jj)+RotorFilletRad2(ii)),0,hpont/2+RotorFilletRad2(ii));
            XpontRadDx(ii) = xS02k(ii);
            YpontRadDx(ii) = hpont/2;
            m = (yS02k(ii))/(xS02k(ii)-x0);
            q = -m*x0;
            [XpontRadBarDx(ii),YpontRadBarDx(ii)] = intersezione_retta_circonferenza(x0,0,-r_all(jj),m,q);

            [xS01k(ii),yS01k(ii)] = intersezione_retta_circonferenza(x0,0,-(r_all(jj+1)-RotorFilletRad1(ii)),0,hpont/2+RotorFilletRad1(ii));
            XpontRadSx(ii) = xS01k(ii);
            YpontRadSx(ii) = hpont/2;
            m = (yS01k(ii))/(xS01k(ii)-x0);
            q = -m*x0;
            [XpontRadBarSx(ii),YpontRadBarSx(ii)] = intersezione_retta_circonferenza(x0,0,-r_all(jj+1),m,q);
        else
            % Se invece i punti 10 e 11 si incrociano (ovvero raccordo troppo grande rispetto)
            % alla larghezza barriera),  disegno solo una linea
            %hpont = 0;
            y10 = pont(ceil(jj/2))/2;
            YpontRadBarDx(ceil(jj/2)) = y10;
            [x_temp,~] = intersezione_retta_circonferenza(x0,0,r_all(jj),0,YpontRadBarDx(ceil(jj/2)));
            XpontRadBarDx(ceil(jj/2)) = 2*x0 - x_temp;
            YpontRadBarSx(ceil(jj/2)) = y10;
            [x_temp,~] = intersezione_retta_circonferenza(x0,0,r_all(jj+1),0,YpontRadBarSx(ceil(jj/2)));
            XpontRadBarSx(ceil(jj/2)) = 2*x0 - x_temp;
            XpontRadDx(ceil(jj/2))  = XpontRadBarDx(ceil(jj/2));
            YpontRadDx(ceil(jj/2))  = y10;
            XpontRadSx(ceil(jj/2))  = XpontRadBarSx(ceil(jj/2));
            YpontRadSx(ceil(jj/2))  = y10;
            RotorFilletRad1(ceil(jj/2)) = 0;
            RotorFilletRad2(ceil(jj/2)) = 0;
        end
    else
        pont(ceil(jj/2)) = 0;
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

% geo.RotorFillet(2,:)=RotorFilletRad;
geo.pontR=pont;
geo.RotorFillet1 = RotorFilletRad1;
geo.RotorFillet2 = RotorFilletRad2;

geo.radial_ribs_split = radial_ribs_split;


