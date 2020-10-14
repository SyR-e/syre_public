%%Calcolo dei punti caratteristici geometria di rotore "Vtype"
%% rev.Gallo 07/05/2018

function [geo,mat,temp]=nodes_rotor_Vtype(geo,mat)

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

p = geo.p;          % Paia poli
nlay = geo.nlay;    % N° layers

dalpha = geo.dalpha;     % Angoli dalpha
pontR=geo.pontR;           % Spessore ponticello radiale 
hc_pu = geo.hc_pu;
dx=geo.dx;

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
XMag5=zeros(1,1); %punti tratto superiore magnete (verso ponticello tangenziale)
YMag5=zeros(1,1);
XMag6=zeros(1,1);
YMag6=zeros(1,1);

XMagpontRadSx=zeros(1,1); %punti tratto inferiore magnete (verso ponticello radiale)
YMagpontRadSx=zeros(1,1);
XMagpontRadDx=zeros(1,1);
YMagpontRadDx=zeros(1,1);

%% DISEGNO DELLA BARRIERA E PONTICELLI RADIALI E TANGENZIALI
%Nota:caso Vtype considero SINGOLA barriera (non è possibile modificare nlay da GUI)
%Definizione punto estremo della barriera
% xpont=(r-pont0)*cos(alpha*pi/180);
% ypont=(r-pont0)*sin(alpha*pi/180);

x0 = r/cos(pi/2/p);
beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,r,x0);
rbeta = (x0 - r * cos(alpha*pi/180))./(cos(beta*pi/180));
[xpont,ypont] = calc_intersezione_cerchi(r-pont0, rbeta, x0);

%% Impostazione inclinazione barriere V-type in  base al parametro "angle" definito nella GUI
%Calcolo punto (Bx0,0) appartenente all'asse della barriera ed all'asse d in base all'angolo di inclinazione angle   
%Escursione massima inclinazione barriera

angle_limit=89*pi/180; %definizione angolo limite di inclinazione barriera [rad]
angle_FEMM=85*pi/180; %angolo massimo inclinazione barriera per definire regione triangolare d'aria sottostante e permettere mesh in FEMM 
if angle > angle_limit %caso limite barriera verticale, setto angolo ad 89° come inclinazione massima limite
    angle=angle_limit; %setto angolo limite a 89°
    disp('#1 Case not allowed, max_angle set at 89°')
end
%Escursione minima inclinazione barriera
if angle < (pi/2/p) %Escursione minima di angle (condizione di parallelismo lato barriera
                    % con lato che delimita regione polo
    angle=pi/2/p;
    disp ('#1 Minimum angle is set to guarantee constant spider')
end

Bx0=xpont-(ypont/tan(angle));

%% Determinazione dello spider e controlli sulla fattibilità della geometria barriera 
% (quanto contenuto dentro function calcHcCheckGeoControlwDx - INIZIO -)

% Spazio di rotore disponibile per l'inserimento della barriera
% Spessori minimi da garantire di ferro (definiti su asse d)
limit_shaft=hfe_min;    %spessore limite ferro rispetto all'albero  
limit_spider=hfe_min; %spessore limite ferro spider
limit_rotor=2*pont0;      %spessore limite ferro su raggio rotore

%la = r - Ar - limit_shaft - limit_rotor - pont0;    % spazio disponibile lungo asse d
la = r - Ar - limit_shaft - limit_rotor;

% Calcolo di max hc e min hc in funzione del parametro alpha
hc_half_min = la/16; %fisso uno altezza minima per metà barriera (1/8 dello spazio disponibile per metà barriera)

% Fisso dei punti limite su asse d (asse di barriera)
Bsh=Ar+limit_shaft;
%Brot=r-pont0-limit_rotor;
Brot=r-limit_rotor;

m_polebound=tan(pi/2/p);
a1=m_polebound; %a1*x+b1*y+c1=0 retta che delimita confine limite barriera verso albero
b1=-1;
c1=-m_polebound*limit_spider;

a2=0;       %a2*x+b2*y+c2=0 retta orizzontale passante per punto estremo barriera
b2=1;
c2=-ypont;
[xBSpider,yBSpider]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
%Distanza a disposizione rispetto a limit_spider
dSpider=calc_distanza_punti([xpont,ypont],[xBSpider,yBSpider]);

%Distanza a disposizione rispetto a limit_shaft
dShaft=calc_distanza_punti([Bsh,0],[Bx0,0]);

xBRot=Brot; %coordinate punto su vincolo limit_rotor avente stessa ordinata del punto estremo
yBRot=ypont;
%Distanza a disposizione rispetto al limite estremo raggio rotore ("lunetta")
dRotor=calc_distanza_punti([xpont,ypont],[xBRot,yBRot]);

hc_half_max=min([dSpider,dShaft,dRotor]); %altezza massima barriera (minima distanza su asse d consentita)
%hc_half_max=la/4;
if hc_half_max < hc_half_min %per casi estremi occorre limitare il valore a quello minimo
    hc_half_max=hc_half_min;
    disp('#2 hc is set to minimum value') 
end

hc = hc_pu * hc_half_max * 2;

%Check controllo valore di hc richiesto da GUI rientrante nei limiti min e max previsti
if hc<2*hc_half_min
    hc=2*hc_half_min;
end
if hc>2*hc_half_max
    hc=2*hc_half_max;
end  

%% Determinazione punti della barriera su asse d nella fasa preliminare
B1k=Bx0-hc/2;
B2k=Bx0+hc/2;

% Inizio dei vari controlli sulle posizioni delle barriere
%% CHECK#1: garantire spider minimo (limit_spider)
m=tan(angle);
m_orto=-1/m;
a1=m; %a1*x+b1*y+c1=0 retta parallela asse barriera passante per B1ktemp
b1=-1;
c1=-m*B1k;
a2=m_orto; %a2*x+b2*y+c2=0 retta perpendicolare asse barriera passante per (xpont,ypont)
b2=-1;
c2=ypont-m_orto*xpont;
[x_temp1,y_temp1]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2); %il check del punto va fatto sulla perpendicolare all'asse della barriera, non sul punto appartenente al cerchio di raggio r-pont0
m_polebound=tan(pi/2/p);
a3=m_polebound; %a3*x+b3*y+c3=0 retta che delimita confine regione polo (polebound)
b3=-1;
c3=0;
[x_temp2,y_temp2]=intersezione_tra_rette(a2,b2,c2,a3,b3,c3);
d=calc_distanza_punti([x_temp1,y_temp1],[x_temp2,y_temp2]); %check spider minimo presente è eseguita 
% su perpendicolare all'asse della barriera stessa 

%Controllo che lato interno della barriera non esca dal confine del polo (alpha molto grandi)
[m1,q1]=retta_abc2mq(a3,b3,c3); 
y_check=m1*x_temp1+q1; %criterio di verifica basato sulla differenza delle ordinate 
if y_check < y_temp1
    disp ('Layer exit to pole bound')
    B1k=B1k+d+limit_spider; %spostamento lato interno barriera verso dx
    hc=(Bx0-B1k)*2; %calcolo del nuovo spessore barriera 
    hc_min=2*hc_half_min;
    if hc < hc_min
        disp('#3 spider is too thin')
    end
elseif d < limit_spider
        
        %Caso in cui lato interno barriera è all'interno della zona di spider minimo
        B1k=B1k+(limit_spider-d); %correzione posizione di B1k su asse d
        disp('#3 spider is too thin')
        hc=(Bx0-B1k)*2; %aggiornamento parametro hc (correggo la geometria per garantire mimimum spider)
        hc_min=2*hc_half_min;
        if hc < hc_min
            disp('#3 hc is less than hc_min to respect minimum spider')
        end
end

%% CHECK#2: garantire minimo spessore ferro su raggio del rotore (limit_rotor)
if (r-B2k) < limit_rotor
    B2k=Brot;
    disp('#4 layer exit from the rotor');
    hc=B2k-B1k; %Aggiornamento parametro hc
    hc_min=2*hc_half_min;
    if hc < hc_min
        disp('#4 hc is less than hc_min to respect limit_rotor')
    end
end

% (quanto contenuto dentro function calcHcCheckGeoControlwDx - FINE -)

%% CHECK 3: Intersezione barriera interna con albero motore (limit_shaft)
% Criterio che uso: cambio la inclinazione della barriera ma mantengo il valore di hc
if (B1k<Bsh)
    B1k=Bsh; %calcolo del B1k limite
    B2k=B1k+hc; %mantengo parametro hc voluto, non lo modifico
    Bx0=(B1k+B2k)/2;  %calcola il nuovo Bx0 limite
    angle_min=atan((ypont/(xpont-Bx0))); %calcolo angle_min che rispetti il nuovo B1k
    disp(['#5 layer barrier cross the shaft, the minimum angle allowed is: ' num2str(angle_min*180/pi)])
    angle=angle_min; %saturo angle ad angle_min (limite inferiore)
end

%Salvataggio parametri che non sono più modificabili nei check successivi
geo.VanglePM=angle; % parametro inclinazione barriere 
geo.hc=hc; % parametro spessore delle barriere 

%% Controllo e correzione della geometria con il grado di libertà dx
dx_old=dx; %salvo il parametro dx introdotto da GUI

%% CHECK#4_Dx Controllo e correzione del dx massimo che posso applicare per evitare l'uscita
% della barriera dalla regione del polo (polebound)
if dx<0 % avoid the exit of the barrier from the pole with dx
    B1ktemp=B1k+dx*hc/2;
    m=tan(angle);
    m_orto=-1/m;
    a1=m; %a1*x+b1*y+c1=0 retta parallela asse barriera passante per B1ktemp
    b1=-1;
    c1=-m*B1ktemp;
    a2=m_orto; %a2*x+b2*y+c2=0 retta perpendicolare asse barriera passante per (xpont,ypont)
    b2=-1;
    c2=ypont-m_orto*xpont;
    [x_temp1,y_temp1]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2); %il check del punto va fatto sulla perpendicolare all'asse della barriera, non sul punto appartenente al cerchio di raggio r-pont0
    
    m_check=tan(pi/2/p);
    a3=m_check; %a3*x+b3*y+c3=0 retta parallela asse barriera ma distante geo.hfemin/2 da polebound (retta che delimita limite area polo, inclinata inclinata di pi/2/p
    b3=-1;
    c3=-m_check*(limit_spider);
    [m1,q1]=retta_abc2mq(a3,b3,c3); 
    y_check=m1*x_temp1+q1;
    
    if y_temp1 > y_check %controllo su ordinata rispetto alla retta limite distante hfe_min/2 da polebound
        [x_temp2,y_temp2]=intersezione_tra_rette(a3,b3,c3,a2,b2,c2); %punto su perpendicolare asse barriera ma appartenente a retta limite barriera interna
        d=calc_distanza_punti([x_temp1,y_temp1],[x_temp2,y_temp2]);
         
        dx1=d*2/hc;
        dx2=abs(dx_old)-dx1;
        dx=-abs(dx2);
       
        disp('#6Dx inner dx modified for prevent exit of the barrier from pole bound')
    end
end

%% CHECK#5_Dx Controllo uscita barriera da rotore con parametro dx
if dx>0
    B2ktemp=B2k+dx*hc/2;
    if (r-B2ktemp)<limit_rotor      
        dx=((Brot-B2k)*2)/hc; %correzione del dx per non violare vincolo limit_rotor
        disp('#7Dx dx modified to prevent exit layer from the rotor')
    end
end

%% CHECK#6_Dx Controllo intersezione barriera con albero motore con parametro dx
if dx<0
    B1ktemp=B1k+dx*hc/2;
    if B1ktemp < Bsh
        dx=((Bsh-B1k)*2)/hc; %correzione dx per non violare vincolo limit_shaft
%         B1k=B1k+dx*hc/2;
%         B2k=B1k+hc; %il valore di Bx0 rimane fisso, non cambia con parametro dx
        disp('#8Dx dx modified to prevent cross the shaft to layer')
    end
end        
    
if dx>0 %ulteriore check per prevenire uscita asse barriera dalla barriera stessa
    B1ktemp=B1k+dx*hc/2;
    if (Bx0-B1ktemp)<geo.pont0 %garantire distanza minima interna asse- lato barriera
        B1ktemp=Bx0-geo.pont0;
        dx1=(B1ktemp-B1k)*2/hc;
        dx=abs(dx1);
        disp('#9Dx dx modified to prevent exit of the axis from the barrier')
    end
elseif dx<0
    B2ktemp=B2k+dx*hc/2;
    if (B2ktemp-Bx0)<geo.pont0 %garantire distanza minima interna asse- lato barriera
        B2ktemp=Bx0+geo.pont0;
        dx2=(B2k-B2ktemp)*2/hc;
        dx=-abs(dx2);
        disp('#10Dx dx modified to prevent exit of the axis from the barrier')
    end
end

%Aggiornamento parametro dx
geo.dx=dx;

%Ri-definizione di B1k e B2k in funzione di dx
B1k = B1k+dx.*hc/2;
B2k = B2k+dx.*hc/2;

%% Calcolo dei punti caratteristici archi di raccordo barriera verso ponticelli tangenziali
% Intersezione circonferenza limite-rette che compongono le barriere 
m2=tan(angle);
q2=-m2.*B2k;
q1=-m2.*B1k;
if angle~=pi/2
    [xTraf2,yTraf2] = intersezione_retta_circonferenza(0,0,r-pont0,m2,q2);
    [xTraf1,yTraf1] = intersezione_retta_circonferenza(0,0,r-pont0,m2,q1);
else
    xTraf1 = B1k;
    yTraf1 = ((r-pont0)^2-xTraf1^2)^0.5;
    xTraf2 = B2k;
    yTraf2 = ((r-pont0)^2-xTraf2^2)^0.5;
end
%Calcolo punti archi di raccordo barriera verso ponticelli tangenziali
%Arco di raccordo sinistro (non soggetto a controlli)
[xt,yt,xc,yc,rc]=tg_cir(B1k,0,xTraf1,yTraf1,xpont,ypont);
xxD1k=xt;
yyD1k=yt;
XcRibTraf1=xc;
YcRibTraf1=yc;

%Arco di raccodo destro (punto di raccordo potrebbe avere ordinata negativa e va limitata a valore nullo)
[xt,yt,xc,yc,rc]=tg_cir(B2k,0,xTraf2,yTraf2,xpont,ypont);
yt_old=yt; %memorizzo ordinata punto di tangenza barriera esterna (servirà per i check successivi sulla barriera esterna)

%CHECK#7:caso in cui ordinata del punto di tangenza yt<0 negativa, barriera esterna composta da solo parte raccordata
if yt<0 %yt<0 (anche se di poco) ma il tratto rettilineo barriera esterna non esiste già più
        %(ho solo parte raccordata per la barriera esterna, alpha<<)
    [xt,yt] = intersezione_retta_circonferenza(xc,yc,rc,0,0);
    B2k = xt; %correzione punto B2k per barriera esterna (punto appartente ad asse d)
end

%Punti caratteristici arco di raccordo sinistro 
xxD2k=xt;
yyD2k=yt;
XcRibTraf2=xc;
YcRibTraf2=yc;
RcRibTraf2=rc;

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
    pontR = F_centrifuga/(sigma_max * l);    % mm
else
    pontR = geo.pont;
end

if (pontR < pont0) % non disegno i ponticelli radiali il cui spessore è minore della tolleranza di lavorazione per gli altri tipi di rotore
    pontR=0;
end
% %% Calcolo dei punti caratteristici archi di raccordo barriera verso ponticelli tangenziali
% % Intersezione circonferenza limite-rette che compongono le barriere 
% m2=tan(angle);
% q2=-m2.*B2k;
% q1=-m2.*B1k;
% [xTraf2,yTraf2] = intersezione_retta_circonferenza(0,0,r-pont0,m2,q2);
% [xTraf1,yTraf1] = intersezione_retta_circonferenza(0,0,r-pont0,m2,q1);
                 
% %Calcolo punti archi di raccordo barriera verso ponticelli tangenziali
% %Arco di raccordo sinistro (non soggetto a controlli)
% [xt,yt,xc,yc,rc]=tg_cir(B1k,0,xTraf1,yTraf1,xpont,ypont);
% xxD1k=xt;
% yyD1k=yt;
% XcRibTraf1=xc;
% YcRibTraf1=yc;

% %Arco di raccodo destro (punto di raccordo potrebbe avere ordinata negativa e va limitata a valore nullo)
% [xt,yt,xc,yc,rc]=tg_cir(B2k,0,xTraf2,yTraf2,xpont,ypont);
% yt_old=yt; %memorizzo ordinata punto di tangenza barriera esterna (servirà per i check successivi sulla barriera esterna)

%CHECK#7:ponticello radiale inserito e yt>0 ancora positivo, barriera esterna composta da solo parte raccordata
%Caso yt>0 risulta ancora positivo ma è al di sotto dello spessore del ponticello radiale
if pontR ~= 0 %check#1: caso con ponticello radiale e ordinata punto di tangenza si trova 
             %al di sotto dello spessore del ponticello pont/2 (alpha<<),
             %MA COMUNQUE yt>0 risulta ancora positivo
    y00=pontR/2;
    if yt<y00 %correzione ordinata punto di tangenza al valori di y00
        [xt,yt] = intersezione_retta_circonferenza(xc,yc,rc,0,y00);
    end         
end    
% %CHECK#8:caso in cui ordinata del punto di tangenza yt<0 negativa, barriera esterna composta da solo parte raccordata
% if yt<0 %yt<0 (anche se di poco) ma il tratto rettilineo barriera esterna non esiste già più
%           %(ho solo parte raccordata per la barriera esterna, alpha<<)
%     if pont==0 %devo distinguere caso con o senza ponticello radiale
%         [xt,yt] = intersezione_retta_circonferenza(xc,yc,rc,0,0);
%         B2k = xt; %correzione punto B2k per barriera esterna (punto appartente ad asse d)
%     else
%         y00=pont/2;
%         [xt,yt] = intersezione_retta_circonferenza(xc,yc,rc,0,y00);
%     end    
% end

%Punti caratteristici definitivi arco di raccordo sinistro 
xxD2k=xt;
yyD2k=yt;
XcRibTraf2=xc;
YcRibTraf2=yc;
RcRibTraf2=rc;

%% Valutazione raccordi barriera ponticello radiale
calc_ribs_rad;

%% Determining  Magnet Area
YcBan = (yyD2k+yyD1k)/2;
XcBan = (xxD2k+xxD1k)/2;
                                                   
Bar_fillfac=geo.BarFillFac;    %barrier filling factor for real magnet

%Disegno area MP- lato superiore rettangolo (verso ponticello tangenziale)
if yyD1k <= yyD2k %devo capire quale punto è più basso per tracciare la retta perpendicolare al lato barriera
                  %per non interferire con uno dei due archi di raccordo
    X5=xxD1k;
    Y5=yyD1k;
    m=tan(angle);
    m_orto=-1/m;
    a1=m_orto;
    b1=-1;
    c1=-m_orto*X5+Y5;
    a2=m;
    b2=-1;
    c2=-m*B2k;
    [X6,Y6]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
    %casi critici (yyD2k molto simile a yyD1k operando con il dx),devo fare un ulteriore 
    %check al criterio precedente 
    %Check definizione lato magnete superiore
    d_check=calc_distanza_punti([B2k,0],[xxD2k,yyD2k]);
    d=calc_distanza_punti([B2k,0],[X6,Y6]);
    if d>d_check
         X6=xxD2k;
         Y6=yyD2k;
         a1=m_orto;
         b1=-1;
         c1=-m_orto*X6+Y6;
         a2=m;
         b2=-1;
         c2=-m*B1k;
         [X5,Y5]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
    end

else
    X6=xxD2k;
    Y6=yyD2k;
    m=tan(angle);
    m_orto=-1/m;
    a1=m_orto;
    b1=-1;
    c1=-m_orto*X6+Y6;
    a2=m;
    b2=-1;
    c2=-m*B1k;
    [X5,Y5]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
    %casi critici (yyD1k molto simile a yyD2k operando con il dx),devo
    %fare un ulteriore check al criterio precedente
    %Check definizione lato magnete superiore
    d_check=calc_distanza_punti([B1k,0],[xxD1k,yyD1k]);
    d=calc_distanza_punti([B1k,0],[X5,Y5]);
    if d>d_check
        X5=xxD1k;
        Y5=yyD1k;
        a1=m_orto;
        b1=-1;
        c1=-m_orto*X5+Y5;
        a2=m;
        b2=-1;
        c2=-m*B2k;
        [X6,Y6]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
        if Y6<0 || Y6<hpont %ulteriore verifica sulla fattibilità area magnete rettangolare
           X6=NaN; %non disegno nulla perchè l'area rettangolare non può essere definita
           Y6=NaN;
        end
       
    end

end

%Salvataggio dei punti calcolati all'interno del vettore apposito
XMag5=X5;
YMag5=Y5;
XMag6=X6;
YMag6=Y6;

%Per la geometria "Vtype" non prevedo alcuna divisione del PM in due parti (vedi caso Circular),
%considero un pezzo unico sfruttando la massima area rettangolare possibile
%all'interno della barriera
if (pontR==0)
    X8=XpontRadBarDx; %punto su barriera lato esterno
    Y8=YpontRadBarDx;
    m=tan(angle);
    m_orto=-1/m;
    a1=m_orto;
    b1=-1;
    c1=-m_orto*X8+Y8;
    a2=m;
    b2=-1;
    c2=-m*B1k;
    [X7,Y7]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
    %Check caso critico su ponticello radiale, definizione lato magnete inferiore (alpha<0.25)
% %     d_check=calc_distanza_punti_altern(B1k,0,xxD1k,yyD1k);
% %     d=calc_distanza_punti_altern(B1k,0,X7,Y7);
% %     if d>d_check %check#1 controllo presenza area magnete rettangolare (valido per alpha <0.25)
% %         X7=NaN; %non disegno nulla perchè area rettangolare magnete non è definibile
% %         Y7=NaN;
% %     end

    %Calcolo atezza magnete da inserire in barriera (anticipo il calcolo per definire in modo 
    %corretto i punti X7,Y7,X8,Y8
    l_check=sqrt((X6-X8)^2+(Y6-Y8)^2); %altezza magnete (se definibile)

    if angle > angle_FEMM && l_check > 0 %check#2 lato inferiore area magnete non deve essere definito, AreaPM è definita
            %AreaPM esiste, lato inferiore area magnete non deve essere disegnato
            X7=B1k; %regione triangolare di aria verso ponticello radiale non definibile
            Y7=0;
            X8=B2k;
            Y8=0;
            disp('Warning: Area Magnet is not rectangular')
    end   
        
    if l_check ==0 % AreaPM non esiste, anche lato superiore area magnete non deve essere disegnato
            XMag5=B1k; %regione triangolare di aria verso ponticello radiale non definibile
            YMag5=0;
            XMag6=B2k;
            YMag6=0;
            X7=B1k; %regione triangolare di aria verso ponticello radiale non definibile
            Y7=0;
            X8=B2k;
            Y8=0;
            disp('Warning: Area Magnet is not defined, only air inside barrier')
    end   
        
                                                                                                                                      
elseif YpontRadBarSx >= YpontRadBarDx %parto a tracciare la perpendicolare dal punto più basso del raccordo con il lato barriera
        X8=XpontRadBarDx;
        Y8=YpontRadBarDx;
        m=tan(angle);
        m_orto=-1/m;
        a1=m_orto;
        b1=-1;
        c1=-m_orto*X8+Y8;
        a2=m;
        b2=-1;
        c2=-m*B1k;
        [X7,Y7]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
    else
        X7=XpontRadBarSx;
        Y7=YpontRadBarSx;
        a1=m_orto;
        b1=-1;
        c1=-m_orto*X7+Y7;
        a2=m;
        b2=-1;
        c2=-m*B2k;
        [X8,Y8]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
end
  
%Salvataggio dei punti calcolati all'interno del vettore apposito
XMagpontRadSx=X7;
YMagpontRadSx=Y7;
XMagpontRadDx=X8;
YMagpontRadDx=Y8;

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
b=sqrt((XMagpontRadSx-XMagpontRadDx)^2+(YMagpontRadSx-YMagpontRadDx)^2); %base area magnete 
l=sqrt((XMag6-XMagpontRadDx)^2+(YMag6-YMagpontRadDx)); %altezza area magnete (da lato esterno);
Area_PM=b*l; %area del magnete

if (angle > angle_FEMM && pontR==0) %porzione di aria triangolare non definibile, area magnete non sarà rettangolare
    
    if Area_PM > 0 %calcolo punto medio solo se Area_PM rettangolare è definibile 
        XmiddleMag_sx=(XMag5+XpontRadBarSx)/2; %punto(su lato sx)retta passante per centro magnete
        YmiddleMag_sx=(YMag5+YpontRadBarSx)/2;

        m=tan(angle);
        m_orto=-1/m;
        a1=m_orto;
        b1=-1;
        c1=-m_orto*XmiddleMag_sx+YmiddleMag_sx;
        a2=m;
        b2=-1;
        c2=-m*B2k;
        [XmiddleMag_dx,YmiddleMag_dx]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
    else
        XmiddleMag_sx=NaN; %punto(su lato sx)retta passante per centro magnete
        YmiddleMag_sx=NaN;
        XmiddleMag_dx=NaN; %punto(su lato dx)retta passante per centro magnete
        YmiddleMag_dx=NaN;
    end
    
else
    if Area_PM > 0 %calcolo area passante per punto medio solo se area PM è positiva
    
        XmiddleMag_sx=(XMag5+XMagpontRadSx)/2; %punto(su lato sx)retta passante per centro magnete
        YmiddleMag_sx=(YMag5+YMagpontRadSx)/2;
        XmiddleMag_dx=(XMagpontRadDx+XMag6)/2; %punto(su lato dx)retta passante per centro magnete
        YmiddleMag_dx=(YMagpontRadDx+YMag6)/2;
    else % non punto medio area rettangolare del magnete
        XmiddleMag_sx=NaN; %punto(su lato sx)retta passante per centro magnete
        YmiddleMag_sx=NaN;
        XmiddleMag_dx=NaN; %punto(su lato dx)retta passante per centro magnete
        YmiddleMag_dx=NaN;
    end
end

if (sum(geo.BarFillFac)==0)
    xmedBar2=(xpont+XMag5+XMag6)/3; %definizione punto medio porzione estrema barriera (parte raccordata superiore verso ponticello tangenziale)
    ymedBar2=(ypont+YMag5+YMag6)/3;
    xair=[xair,xmedBar2];
    yair=[yair,ymedBar2];
    xmag=[xmag,cosd(0)];
    ymag=[ymag,sind(0)];
    Br = [Br mat.LayerMag.Br];
    
else
    if Area_PM >0
        % magnetization direction  BarfillFac~=0
        [a1,b1,c1]=retta_per_2pti(XmiddleMag_sx,YmiddleMag_sx,XmiddleMag_dx,YmiddleMag_dx);
        [a2,b2,c2]=retta_per_2pti(Bx0,0,XcBan,YcBan);

        [xPM,yPM]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2); %punto centrale posizione magnete

        angle_PM=atand(-a1/b1); %direzione vettore induzione residua magnete (orientamento vettore Br uscente da rotore, non ancora ruotato come da impostazione SPM)

        %Definizione punti medi e proprietà porzioni di barriera
        xmedBar1=xPM; %punto medio area magnete PM
        ymedBar1=yPM;

        if strcmp (mat.LayerMag.MatName,'Bonded-Magnet') || strcmp(mat.LayerMag.MatName,'Air')
            
            xair=[xair,xmedBar1];  %vettore punti medi aree di aria o magnete fittizio
            yair=[yair,ymedBar1];

            Br = [Br mat.LayerMag.Br]; %valore di induzione residua dell'aria

            xmagair=[xmagair,cosd(0)]; %orientamento vettore di induzione residua magnete fittizio o aria 
            ymagair=[ymagair,sind(0)]; 

        else
            
            xc=[xc,xmedBar1]; %aggiornamento vettore punti medi porzioni barriere di magnete
            yc=[yc,ymedBar1];

            Br = [Br mat.LayerMag.Br]; %valore di induzione residua da assegnare alla porzione di magnete

            xmag=[xmag,-cosd(angle_PM)]; %definzione orientamento vettore di induzione residua Br del magnete (criterio usato: perpendicolare al lato interno barriera)
            ymag=[ymag,-sind(angle_PM)]; % vettore Br definito come entrante nel rotore (come da impostazione SPM)
        end
    end
    
    %Definisco porzione di aria barriera verso ponticello tangenziale
    xmedBar2=(xpont+XMag5+XMag6)/3; %definizione punto medio porzione estrema barriera (parte raccordata superiore verso ponticello tangenziale)
    ymedBar2=(ypont+YMag5+YMag6)/3;
        
    if pontR==0 %distinzione tra caso con e senza ponticelli radiali per inviduare le regioni
                %da assegnare i labels                
        if angle > angle_FEMM %porzione di aria sotto il magnete non è definibile, per cui non la disegno
            
            if Area_PM > 0 %area magnete ben definita, non disegno solo triangolino d'aria sottostante         
                xair=[xair,xmedBar2];   %vettore punti medi aree di aria o magnete fittizio
                yair=[yair,ymedBar2];

                xmagair=[xmagair,cosd(0)]; %orientamento vettore di induzione residua dell'aria 
                ymagair=[ymagair,sind(0)];

                Br = [Br mat.LayerMag.Br]; %valore di induzione residua dell'aria
                
            else %area rettangolare magnete non definibile, solo triangolino d'aria (sul lato esterno non ho tratto lineare)
                xmedBar3=(B1k+B2k+XMag5)/3; 
                ymedBar3=YMag5/3;
                
                xair=[xair,xmedBar2,xmedBar3]; %aggiornamento vettore punti medi aree di aria o magnete fittizio
                yair=[yair,ymedBar2,ymedBar3];

                Br = [Br mat.LayerMag.Br mat.LayerMag.Br]; %valore di induzione residua dell'aria

                xmagair=[xmagair,cosd(0),cosd(0)]; %orientamento vettore di induzione residua dell'aria 
                ymagair=[ymagair,sind(0),sind(0)];
            end
                
        else %Caso in cui il triangolino di aria al di sotto del magnete è definibile     
            
            %Caso senza ponticello radiale:la regione d'aria al di sotto del magnete è sempre triangolare       
            xmedBar3=(B1k+B2k+XMagpontRadSx)/3; 
            ymedBar3=YMagpontRadSx/3;
                
            xair=[xair,xmedBar2,xmedBar3]; %aggiornamento vettore punti medi aree di magnete
            yair=[yair,ymedBar2,ymedBar3];

            Br = [Br mat.LayerMag.Br mat.LayerMag.Br]; %valore di induzione residua dell'aria

            xmagair=[xmagair,cosd(0),cosd(0)]; %orientamento vettore di induzione residua dell'aria 
            ymagair=[ymagair,sind(0),sind(0)];
                   
         end 
        
    else
    %Caso con ponticelli radiali: occorre distinguere due casi, regione
    %triangolare e regione approssimabile come rettangolare
        if yt_old<0 || (yt_old>0 && yt_old<hpont) %regione di aria sottostante è triangolare
        xmedBar3=(x_temp1+XMagpontRadSx+xxD2k)/3;
        ymedBar3=(y_temp1+YMagpontRadSx+yyD2k)/3;
        else
            XmiddleAir2_sx=(XMagpontRadSx+x_temp1)/2;
            YmiddleAir2_sx=(YMagpontRadSx+y_temp1)/2;

            XmiddleAir2_dx=(XMagpontRadDx+x_temp2)/2;
            YmiddleAir2_dx=(YMagpontRadDx+y_temp2)/2;

            [a1,b1,c1]=retta_per_2pti(XmiddleAir2_sx,YmiddleAir2_sx,XmiddleAir2_dx,YmiddleAir2_dx);
            [a2,b2,c2]=retta_per_2pti(Bx0,0,XcBan,YcBan);

            [xmedBar3,ymedBar3]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
        end
        
        xair=[xair,xmedBar2,xmedBar3];
        yair=[yair,ymedBar2,ymedBar3];
        xmagair=[xmagair,cosd(0),cosd(0)];
        ymagair=[ymagair,sind(0),sind(0)];
        Br = [Br mat.LayerMag.Br mat.LayerMag.Br];       
    end
end

if sum(geo.BarFillFac)==0
    mat.LayerMag.Br = [mat.LayerMag.Br mat.LayerMag.Br];    % replicates Br for correct block assignation
else
    % ho 2 segmenti per ogni mezzo layer
    BrTemp=[];
    for ii=1:length(mat.LayerMag.Br)
        BrTemp=[BrTemp,mat.LayerMag.Br(ii) mat.LayerMag.Br(ii)];
    end
    mat.LayerMag.Br = BrTemp;
end

%Salvataggio dei punti che individuano le varie aree all'interno della
%barriera (aria o magnete commerciale precedentemente scelto)
temp.xc=xc;
temp.yc=yc;
temp.xmag=xmag;
temp.ymag=ymag;
temp.zmag=zeros(1,length(xmag));

temp.xair = xair;
temp.yair = yair;
temp.xmagair= xmagair;
temp.ymagair= ymagair;

%% Salvataggio dei dati finali:
geo.pont = pontR;
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
temp.XcRibTraf1=XcRibTraf1;
temp.YcRibTraf1=YcRibTraf1;
temp.XcRibTraf2=XcRibTraf2;
temp.YcRibTraf2=YcRibTraf2;
temp.RcRibTraf2=RcRibTraf2;

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
temp.XMag5=XMag5;
temp.YMag5=YMag5;
temp.XMag6=XMag6;
temp.YMag6=YMag6;

temp.XMagpontRadSx=XMagpontRadSx;
temp.YMagpontRadSx=YMagpontRadSx;
temp.XMagpontRadDx=XMagpontRadDx;
temp.YMagpontRadDx=YMagpontRadDx;

%% Center points of the PM segments
if (sum(geo.BarFillFac)~=0 && Area_PM >0)
    temp.xPM = xPM;
    temp.yPM = yPM;
end

end