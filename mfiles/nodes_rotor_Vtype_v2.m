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

function [geo,mat,temp]=nodes_rotor_Vtype_v2(geo,mat)

%% Inizializzazione dei dati in ingresso
angle=geo.VanglePM;    % Inclinazione barriera [radianti]
r = geo.r;          % Raggio del rotore al traferro
R = geo.R;          % Raggio esterno totale
lt = geo.lt;        % Lunghezza denti

x0 = geo.x0;        % Centro fittizio
rshaft = geo.Ar;    % Raggio albero
Ar=geo.Ar;
l = geo.l;          % Lunghezza pacco lamierini
g = geo.g;          % Spessore traferro
pont0 = geo.pont0;  % Spessore ponticelli al traferro (minima tolleranza meccanica di lavorazione)
pontT = geo.pontT(1);  % Airgap ribs thickness [mm]
p = geo.p;          % Paia poli
nlay = geo.nlay;    % N° layers

dalpha = geo.dalpha;     % Angoli dalpha
%pont=geo.pont(1:nlay);           % Spessore ponticello radiale 
hc_pu = geo.hc_pu;
%dx=geo.dx;

racc_pont = geo.racc_pont;      % racc_pont=1*pont0 <- per i ponticelli radiali.

nmax = geo.nmax;                % Velocità max (rpm) per la valutazione della sollecitazione centrifuga più gravosa (-> vedi ponticelli radiali)
hfe_min=geo.hfe_min;            % Spessore minimo del ferro per garantire passaggio delle linee di flusso guida
sigma_max = mat.Rotor.sigma_max;    % Snervamento materiale [MPa]
rhoFE = mat.Rotor.kgm3;             % Densità del ferro di rotore [kg/m3]
rhoPM = mat.LayerMag.kgm3;          % Densità magneti [kg/m3]

% Valutazione alpha (nlay=1)
alpha= dalpha;      % angolo di riferimento posizione asse barriera 

% Inizializzazione punti caratteristici della barriera - parte superiore
% barriera (ragiono su mezzo layer)

XcRibTraf1=zeros(1,1); %punti riguardanti archi di raccordo verso ponticelli tangenziali
XcRibTraf2=zeros(1,1);
YcRibTraf1=zeros(1,1);
YcRibTraf2=zeros(1,1);
xxD1k=zeros(1,1);
yyD1k=zeros(1,1);
xxD2k=zeros(1,1);
yyD2k=zeros(1,1);

xTraf1=zeros(1,1);
yTraf1=zeros(1,1);
xTraf2=zeros(1,1);
yTraf2=zeros(1,1);

% Inizializzazione punti caratteristici della barriera - parte inferiore
% barriera verso asse d (ragiono su mezzo layer)
Bx0=zeros(1,1);

% Inizializzazione punti individuazione area magnete (se questo è presente)
xPMC1t=zeros(1,1); %punti tratto superiore magnete (verso ponticello tangenziale)
yPMC1t=zeros(1,1);
xPMC2t=zeros(1,1);
yPMC2t=zeros(1,1);

xPMC1b=zeros(1,1); %punti tratto inferiore magnete (verso ponticello radiale)
yPMC1b=zeros(1,1);
xPMC2b=zeros(1,1);
yPMC2b=zeros(1,1);

%% DISEGNO DELLA BARRIERA E PONTICELLI RADIALI E TANGENZIALI
%Nota:caso Vtype considero SINGOLA barriera (non è possibile modificare nlay da GUI)
%Definizione punto estremo della barriera
xpont=(r-pontT).*cos(alpha*pi/180);
ypont=(r-pontT).*sin(alpha*pi/180);

% x0 = r/cos(pi/2/p);
% beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,r,x0);
% rbeta = (x0 - r * cos(alpha*pi/180))./(cos(beta*pi/180));
% [xpont,ypont] = calc_intersezione_cerchi(r-pont0, rbeta, x0);

%% Impostazione inclinazione barriere V-type in  base al parametro "angle" definito nella GUI
%Calcolo punto (Bx0,0) appartenente all'asse della barriera ed all'asse d in base all'angolo di inclinazione angle   
%Escursione massima inclinazione barriera

angle_limit=90*pi/180; %definizione angolo limite di inclinazione barriera [rad]
angle_FEMM=85*pi/180; %angolo massimo inclinazione barriera per definire regione triangolare d'aria sottostante e permettere mesh in FEMM 
if angle > angle_limit %caso limite barriera verticale, setto angolo ad 89° come inclinazione massima limite
    angle=angle_limit; %setto angolo limite a 89°
    disp(['#1 Case not allowed, max_angle set at ' int2str(angle_limit) ' deg'])
end
%Escursione minima inclinazione barriera
if angle < (pi/2/p) %Escursione minima di angle (condizione di parallelismo lato barriera
                    % con lato che delimita regione polo
    angle=pi/2/p;
    disp ('#1 Minimum angle is set to guarantee constant spider')
end



%% Determinazione dello spider e controlli sulla fattibilità della geometria barriera 
% (quanto contenuto dentro function calcHcCheckGeoControlwDx - INIZIO -)

% Spazio di rotore disponibile per l'inserimento della barriera
% Spessori minimi da garantire di ferro (definiti su asse d)
limit_shaft=hfe_min;    %spessore limite ferro rispetto all'albero  
limit_spider=hfe_min; %spessore limite ferro spider
limit_rotor=pont0;      %spessore limite ferro su raggio rotore

cos_x0 = cos(pi/2/p);
sin_x0 = sin(pi/2/p);
ArLim = geo.r*(1/cos_x0-sqrt(((1-cos_x0^2)/cos_x0)^2+sin_x0^2));    % shaft radius according to the pole span and x0

%la = r - Ar - limit_shaft - limit_rotor - pont0;    % spazio disponibile lungo asse d
% la = r - Ar - limit_shaft - limit_rotor;
la = r - ArLim - limit_shaft - limit_rotor;

% Calcolo di max hc e min hc in funzione del parametro alpha
% hc_half_min = 2*pont0; % fisso uno altezza minima per metà barriera lungo l'asse d (1/8 dello spazio disponibile per metà barriera)
% hc_half_max=la/4;       % altezza massima della barriera lungo l'asse d pari a metà dello spazio disponibile

%hc = hc_pu * hc_half_max * 2;

% % SF: normalizzazione cambiata dicembre 2018
hc_half_min = pont0;
hc_half_max = la/4;
if hc_pu<1
    hc = (hc_half_min+hc_pu*(hc_half_max-hc_half_min))*2;
    warning('hc_pu correct with the new method!!! Save the machine!!!')
    hc_pu=hc/g;
end
% % con questa normalizzazione hc_pu corretto è compreso sempre tra 0 e 1
% % (limiti di alpha a parte)

hc = hc_pu*g;

%Check controllo valore di hc richiesto da GUI rientrante nei limiti min e max previsti
if hc<2*hc_half_min
    hc=2*hc_half_min;
    disp('#1: min hc set')
end

if hc>2*hc_half_max
    hc=2*hc_half_max;
    disp('#1: max hc set')
end

hc_pu=hc/g;


%% Check #1: hc <vs> spider
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
yS = xS.*m2+q2;
hcMaxSpider = 2*((xcbar-xS).^2+(ycbar-yS).^2).^0.5;

if hc>hcMaxSpider
    disp('#1 thin spider, hc limit')
    hc=hcMaxSpider;
    %hc_pu=hc/(2*hc_half_max);
    hc_pu = (hc/2-hc_half_min)/(hc_half_max-hc_half_min);
    xcbar = (r-pontT-hc/2).*cos(alpha_pont);
    ycbar = (r-pontT-hc/2).*sin(alpha_pont);
end

% update hc in geo
geo.hc = hc;
geo.hc_pu = hc_pu;

%% Determinazione punti della barriera su asse d nella fasa preliminare
Bx0=xcbar-(ycbar/tan(angle));
B1k=Bx0-hc/2/sin(angle);
B2k=Bx0+hc/2/sin(angle);

%% Check #2: hc <vs> shaft
intShaft=B1k-Ar-hfe_min;
if intShaft<0
    disp('#2 barrier cross shaft, angle limit')
    angle = calc_Vtype_angle_lim(geo);
    Bx0   = xcbar-(ycbar/tan(angle));
    B1k   = Bx0-hc/2/sin(angle);
    B2k   = Bx0+hc/2/sin(angle);
end

% update hc in geo
geo.VanglePM = angle;

%% Calcolo dei punti caratteristici archi di raccordo barriera verso ponticelli tangenziali
% retta 1: asse della barriera
% retta 2: perpendicolare a retta 1 e passante per il centro dell'arco terminale

m1 = (ycbar)/(xcbar-Bx0);
q1 = -m1*Bx0;
m2 = -1/m1;
q2 = ycbar-m2*xcbar;
xxD1k = xcbar-abs(hc/2*cos(atan(m2)));
yyD1k = xxD1k*m2+q2;
xxD2k = xcbar+abs(hc/2*cos(atan(m2)));
yyD2k = xxD2k*m2+q2;

%% Valutazione ponticelli radiali: (criterio usato simile al caso Seg)- rev.Gallo 07/05/2018
% Nota: considero anche la massa appesa dei magneti da rev.327 (Gallo)                
% Calcolo spessore ponticello radiale, distinguendo caso inserimento manuale o automatico in base alla velocità di rotazione                 
Afe=zeros(1,1);
tmp_bary=zeros(1,2);

rTemp=abs((x0-xpont)+1j*ypont);

% Calcolo area Ferro
[xrot_traf,yrot_traf]=calc_intersezione_cerchi(r, rTemp, x0);
Dx=(r-xrot_traf(1))/5;  %per avere almeno 5 divisioni;
xcir_plot=[r:-Dx:r*cos(pi/2/p)];
ycir_plot=sqrt(r^2-xcir_plot.^2);
VectCir=find(xcir_plot>=xrot_traf(1));
x_ext_rot=xcir_plot(VectCir);
y_ext_rot=ycir_plot(VectCir);

X=[B2k, xxD2k,xpont,fliplr(x_ext_rot)];
Y=[0,yyD2k,ypont,fliplr(y_ext_rot)];
Afe=polyarea(X,Y); %area del poligono identificato dai punti di X e Y
tmp_bary=centroid(X',Y'); %centro di massa serie di punti
clear X Y;

rG=tmp_bary(1);     % baricentro della regione di ferro sorretta dal ponticello a I

M_Fe = 2*Afe*l * 1e-9 * rhoFE ;   % massa ferro appeso al ponticello

%Calcolo area magnete (barriera a "V")
X_Ibarr=[xxD1k xxD2k B2k B1k];
Y_Ibarr=[yyD1k yyD2k 0 0];

areaV=2*polyarea(X_Ibarr,Y_Ibarr);  % area della barriera a V in cui inserire i magneti (escludo i due semicerchi all'estremità della barriera e gli eventuali ponticelli)
A_PM=areaV;
M_PM=A_PM*rhoPM*1e-9*2*l; % massa magnete appesa al ponticello

%Calcolo spessore ponticelli radiali (tengo conto area ferro+area magnete)
F_centrifuga = (M_Fe+M_PM) .* rG/1000 *  (nmax * pi/30)^2;
        
if geo.radial_ribs_eval == 0
    pont = F_centrifuga/(sigma_max * l);    % mm
else
    pont = geo.pontR;
end

if (pont < pont0) % non disegno i ponticelli radiali il cui spessore è minore della tolleranza di lavorazione
    pont=0;
end


%% Valutazione raccordi barriera ponticello radiale
yt_old=yyD2k;
calc_ribs_rad;

%% Determining  Magnet Area
YcBan = (yyD2k+yyD1k)/2;
XcBan = (xxD2k+xxD1k)/2;

%Salvataggio dei punti calcolati all'interno del vettore apposito
xPMC1t=xxD1k;
yPMC1t=yyD1k;
xPMC2t=xxD2k;
yPMC2t=yyD2k;

%Per la geometria "Vtype" non prevedo alcuna divisione del PM in due parti (vedi caso Circular),
%considero un pezzo unico sfruttando la massima area rettangolare possibile
%all'interno della barriera
X8=XpontRadBarDx;
Y8=YpontRadBarDx;
X7=X8-hc*sin(angle);
Y7=Y8+hc*cos(angle);
if ((angle > angle_FEMM)&&(pont==0)) %check#2 su angle_FEMM, se vero, il magnete non è rettangolare e il bordo inferiore è l'asse d (mesh FEMM)
    X7=B1k;
    Y7=0;
    X8=B2k;
    Y8=0;
    %disp('Warning: PM shape is not rectangular')
end 
  
%Salvataggio dei punti calcolati all'interno del vettore apposito
xPMC1b=X7;
yPMC1b=Y7;
xPMC2b=X8;
yPMC2b=Y8;

%% Determinazione dei differenti segmenti magnete,punto centrale (baricentro geometrico)
% e direzione di magnetizzazione vettore induzione residua Br

%Inizializzazione variabili necessarie ad identificare le varie aree all'interno della barriera
xmag=[];
ymag=[];
xair=[];
yair=[];
Br = [];
xc=[];
yc=[];
xmagair=[];
ymagair=[];

%Determinazione punto medio area magnete rettangolare (se tale regione esiste)
%CHECK #1: calcolo area magnete rettangolare (per capire se è definibile o no)
b=sqrt((xPMC1b(1)-xPMC2b(1))^2+(yPMC1b(1)-yPMC2b(1))^2); %base area magnete 
l=sqrt((xPMC2t(1)-xPMC2b(1))^2+(yPMC2t(1)-yPMC2b(1))); %altezza area magnete (da lato esterno);
Area_PM=b*l; %area del magnete
% punti centrali magnete
xc=(xxD1k+xxD2k+X7+X8)/4;
yc=(yyD1k+yyD2k+Y7+Y8)/4;
xPM=xc;
yPM=yc;
angle_PM=angle+pi/2;
xmag=cos(angle_PM);
ymag=sin(angle_PM);
% punti centrali semicerchi di aria
xair=(xpont+xPMC1t+xPMC2t)/3; %definizione punto medio porzione estrema barriera (parte raccordata superiore verso ponticello tangenziale)
yair=(ypont+yPMC1t+yPMC2t)/3;
if pont>0
    xtemp=(XpontRadSx+XpontRadBarSx+XpontRadDx+XpontRadBarDx)/4;
    ytemp=(YpontRadSx+YpontRadBarSx+YpontRadDx+YpontRadBarDx)/4;
    xair=[xair,xtemp];
    yair=[yair,ytemp];
elseif ((angle<=angle_FEMM)&&(angle~=90)) % magnete rettangolare, devo aggiungere porzione di aria
    xtemp=(B1k+B2k+xPMC1b)/3; 
    ytemp=(0+0+yPMC1b)/3;
    xair=[xair,xtemp];
    yair=[yair,ytemp];
end
xmagair=zeros(size(xair));
ymagair=xmagair;




%Salvataggio dei punti che individuano le varie aree all'interno della
%barriera (aria o magnete commerciale precedentemente scelto)
temp.xcbar=xc;
temp.ycbar=yc;
temp.xmag=xmag;
temp.ymag=ymag;
temp.zmag=zeros(size(xmag));

temp.xair = xair;
temp.yair = yair;
temp.xmagair= xmagair;
temp.ymagair= ymagair;

%% Salvataggio dei dati finali:
geo.pontR = pont;
geo.hc=hc;
geo.VanglePM=angle;
geo.B1k=B1k;
geo.B2k=B2k;
geo.Bx0=Bx0;

temp.B1k=B1k;
temp.B2k=B2k;
temp.Bx0=Bx0;
temp.xpont=xpont;
temp.ypont=ypont;
temp.xxD1k=xxD1k;
temp.yyD1k=yyD1k;
temp.xxD2k=xxD2k;
temp.yyD2k=yyD2k;
temp.XcRibTraf1=xcbar;
temp.YcRibTraf1=ycbar;
temp.XcRibTraf2=xcbar;
temp.YcRibTraf2=ycbar;
temp.RcRibTraf2=hc/2;
temp.xc=xc;
temp.yc=yc;

%% Points for radial ribs
temp.XpontRadDx=XpontRadDx;
temp.YpontRadDx=YpontRadDx;
temp.XpontRadSx=XpontRadSx;
temp.YpontRadSx=YpontRadSx;
temp.XpontRadBarDx=XpontRadBarDx;
temp.XpontRadBarSx=XpontRadBarSx;
temp.YpontRadBarDx=YpontRadBarDx;
temp.YpontRadBarSx=YpontRadBarSx;

%% Centro e raccordi ponticelli radiali
temp.XcRaccpontRadSx=XcRaccpontRadSx;
temp.YcRaccpontRadSx=YcRaccpontRadSx;
temp.RcRaccpontRadSx=RcRaccpontRadSx;
temp.XcRaccpontRadDx=XcRaccpontRadDx;
temp.YcRaccpontRadDx=YcRaccpontRadDx;
temp.RcRaccpontRadDx=RcRaccpontRadDx;

%% Points for magnet segmentation (BarFillFact~=0)
temp.xPMC1t = xPMC1t;
temp.yPMC1t = yPMC1t;
temp.xPMC2t = xPMC2t;
temp.yPMC2t = yPMC2t;
temp.xPMC1b = xPMC1b;
temp.yPMC1b = yPMC1b;
temp.xPMC2b = xPMC2b;
temp.yPMC2b = yPMC2b;

%% Center points of the PM segments
if (Area_PM >0)
    temp.xPM = xPM;
    temp.yPM = yPM;
end

mat.LayerMag.Br = [mat.LayerMag.Br mat.LayerMag.Br];

end