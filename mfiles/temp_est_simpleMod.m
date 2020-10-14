% Soluzione rete semplificata, riferimento formule articolo Simplified
% [Boglietti et all,Thermal model for variable speed self cooled industrial
% Induction Motor]
%

function [Temp]=temp_est_simpleMod(geo,per)
%% Thermal material initialization
% calore specifico carcassa (frame) [J/g/K]
Materiali.Carcassa.Cth=1005.6;
% Resistenza termica carcassa
Materiali.Carcassa.Rth=0.025;
% calore specifico flangia (Endcap) [J/kg/K]
Materiali.Carcassa.Cth=1005.6;
% densità di carcassa
Materiali.Carcassa.d=2700;
% Conducibilita del ferro di statore in laminazione assiale [W/m K] + MANUALE
Materiali.FerroSta.Rth=25;
% Conducibilità termica radiale del ferro di statore [W/m K]
Materiali.FerroSta.Rth=25;
% Densità del ferro
Materiali.FerroSta.d=7800;
% calore specifico del ferro di statore [J/kg/K]
Materiali.FerroSta.Cth=490;
Materiali.Conduttore.Rth=387;
% calore specifico del rame [J/kg/K]
Materiali.Conduttore.Cth=393;
% densità del rame [kg/m3]
Materiali.Conduttore.d=8500;
% coefficiente conduttività film di traferro (statico o dinamico) [W/m^2 K]
Materiali.Isolante.Rth=0.033;
% calore specifico aria
Materiali.Isolante.Cth=1014;
% conducibilità termica dell'acciaio dell'albero [W/m K]
Materiali.Albero.Rth=1;
% densità ferro albero
Materiali.Albero.d=7800;
% capacità termica del acciaio albero [J/kg/K]
Materiali.Albero.Cth=1;
% densità aria
Materiali.Isolante.d=1.2;
%
% ManTherm.R1     = 0.0328; 
ManTherm.hc     = 1*10^3; % frame-core contact coefficient (1991-Mellor-LPTM...)
%ManTherm.hc     = 300;

% R1 = ManTherm.R1; 
hc = ManTherm.hc;

%% Data Initialization
Pjs=per.Loss;
tetaFrame = per.temphous;
% max speed or rated speed
omega=geo.nmax;     %rpm
% geometrical data
geo.ly = (geo.R - (geo.r + geo.g + geo.lt))/geo.r;

n3ph = geo.win.n3phase;
p    = geo.p;
q    = geo.q;
Kcu  = geo.win.kcu;
g    = geo.g*10^-3;
R    = geo.R*10^-3;
L    = geo.l*10^-3;
wt   = geo.wt*10^-3;
lt   = geo.lt*10^-3;
ly   = geo.ly*10^-3;
x    = geo.r/geo.R;
r    = geo.r*1e-3;
ps   = 2*pi*x*R/(6*p*q*n3ph);    % slot pitch (m)

% Thermal properties
LambdaFE=Materiali.FerroSta.Rth;
lamdaEQ=0.3;
LambdaNomex=0.14;

%% %%%%%%%%%%%%%%%%%
%%  End winding Calc
%% %%%%%%%%%%%%%%%%%
% if q>=1
%     wt=2*pi*kt*x*b./(6*p*q)*R;
%     ly=b*x*R./(p);
%     lt=R-ly-x*R-g;
%     ltestata=(2*lt+kracc*(0.5*pi*R*(1-ly+x)./p));
%
% else
%     wt=2*pi*kt*x*b./(6*p*q)*R;
%     ly=b*x*R./(p);
%     lt=R-ly-x*R;
%     kt=wt*6*p*q/(2*pi*x*b*R);
%     ltestata=0.5*(wt+pi*(x*R+lt/2).*sin(pi./(6*p*q)));
%
% end

% Kend=(ltestata+L)/L;

%% %%%%%%
%% end winding end
%%  %%%%%

%% %%%%%%
%% Thermal resistance calc
%%  %%%%%
% R2 core to frame contact resistace
R2 = 1/(2*pi*hc*L*R);
R2 = 1/(pi*hc*L*R);
%
% lN=0.0002;
% ws=ps-wt;
% Rse=x*R+lt;
% Aslot=(pi*R^2*(1-b*x/p).*(1+x-2*b*kt*x))/(6*p*q)
if isfield(geo,'Aslot')
    Aslot=geo.Aslot*1e-6;%/(6*pi*q);
else
    Aslot=2*pi/(6*p*q*n3ph)*lt*(R-lt/2-ly)-wt*lt;
end
Acu=Aslot*Kcu;
Rsa=r+g;
Rsb=Rsa+lt;
wsEq=0.5*(Rsa+Rsb)*pi/3/p/q/n3ph-wt;
hsEq=Aslot./wsEq;
% hsEq=lt;
% wsEq=Aslot/lt;

% wsp=wsEq.*(1-Kcu);

Ps=2*hsEq+wsEq;

teq=((1-Kcu)*Aslot./Ps);
% Rcu_ir=teq./lamdaEQ./(Ps*L)/(6*p*q*n3ph);
%
n_element=2;
for kk=1:n_element
    eval(['Rcu',num2str(kk),'=0.5*teq./lamdaEQ./(2*hsEq/',num2str(n_element),'*L)/(6*p*q*n3ph);']);
end
eval(['Rcu',num2str(n_element+1),'=0.5*teq./lamdaEQ./(wsEq*L)/(6*p*q*n3ph);']);
%
% Rcu11=teq./lamdaEQ./(hsEq/2*L)/(6*p*q);
% Rcu22=Rcu11;
% Rcu33=teq./lamdaEQ./(wsEq*L)/(6*p*q);
%
w=(wsEq-2*teq)/3;
h1=hsEq-teq-w;
Acu1=w*h1;
Ps1=2*hsEq-w;
Rcu01=0.5*teq/(lamdaEQ*Ps1*L)/(6*p*q*n3ph);
% yoke
Rsy=1/(2*pi*LambdaFE*L)*log(R./Rsb);
% Rsy = 1/(2*pi*LambdaFE*L)*R/(r+g+lt);

% teeth
pir=(6*p*q*n3ph)*wt.*lt./((6*p*q*n3ph)*wt.*lt+(6*p*q*n3ph)*Aslot);
% Rst=1./(2*pi*LambdaFE*L*pir).*log(Rsb./Rsa);

% rm=(Rsa+Rsb)/2;
% fp=2*pi/(6*p*q*n3ph);
% fe=wt./Rsa;

for kk=1:n_element
    eval(['Rst',num2str(kk),'=1./(2*pi*LambdaFE*L*pir).*log((Rsa+lt*',num2str(kk/n_element),')./(Rsa+lt*',num2str((kk-1)/n_element),'));']);
%     eval(['Rst',num2str(kk),'=1./(2*pi*LambdaFE*L*pir).*((Rsa+lt*',num2str(kk/n_element),')./(Rsa+lt*',num2str((kk-1)/n_element),'));']);
end

% Rst11=1./(2*pi*LambdaFE*L*pir).*log((Rsa+lt/2)./Rsa);
% Rst22=1./(2*pi*LambdaFE*L*pir).*log((Rsa+lt)./(Rsa+lt/2));

% Rst1 = 1/(2*pi*LambdaFE*pir*L)*(r+g+lt/2)/(r+g);
% Rst2 = 1/(2*pi*LambdaFE*pir*L)*(R)/(r+g+lt/2);


%
Rew_ec=1./(2*pi*0.2*lt).*log(R./(R-0.5*ly));
% Rgap
eta=(Rsa-g)./Rsa;
RE=Rsa*omega*pi/30*g/(20*10^-6);
lamGAP=0.0019*eta.^-2.9084+RE.^(0.4614*log(0.3361*eta));
% lamGAP=0.13;
Rgap=1./(2*pi*L*lamGAP).*log(Rsa./(Rsa-g));
%
Rsig=R2;
% Rew,ia
Aew=1.1*lt*2*pi.*Rsa;
v=omega*pi/30*x*R*0.5;
hew=15.5*(0.29*v+1);
Rew_ia=1./(Aew.*hew);
% Ria,ec
Aec=2*pi*(R+0.0106)^2;
Ria_ec=1./(Aec.*hew);
% Rs,ag
Aist=2*pi*Rsa*L;
Nu=2;
% lambdaAIR=Materiali.Isolante.Rth;
lambdaAIR=0.026;

%% calcolo dei numeri per il traferro
% calcolo dei numeri di Taylor e Prandtl TODO 
% [Nut,NuEw] = NusseltNumb(0);

Nta = 1;
Npr = 1;

if (0)    
    if (Nta <= 41)
        Nut = 2.2;
    else
        Nut = 0.23 * (Nta^0.63) * (Npr^0.27);
    end
else
    Nut = 2.0;
end

%% calcolo del numero di Nusselt per le connessionni frontali
% numero di Reynolds e di Rayleigh TODO
% Nre = 1;
% Nra = 1;
% NuEw = 0.83*(Nre*Nra)^0.2;

hag = Nut*lambdaAIR/g;
Rsag=1./(Aist.*hag);
% Rrag
Rrag=1./(2*pi*x.*R*L*hag);
% Rshaft
% Rshaft=0.0645;
% rotor +shaft
Rrshf=1./(2*pi*LambdaFE*L)*log(1.5*x*R./(x*R-ly))+0.25*0.5*L./(LambdaFE*pi*0.08^2)+0.5*0.5*1.3*L./(LambdaFE*pi*0.08^2);
% %%%
%     Req1=Rst+Rsy+Rcu_ir;
%     Rtmp=(Rcu11+Rst11).*Rcu22./(Rcu11+Rst11+Rcu22)+Rst22;
%     
%     for kk=1:n_element
%         if kk==1
%             Rtemp=((Rcu1+Rst1).*Rcu2)./(Rcu1+Rst1+Rcu2);
%         else
%             eval(['Rtemp=((Rtemp+Rst',num2str(kk),').*Rcu',num2str(kk+1),')./(Rtemp+Rst',num2str(kk),'+Rcu',num2str(kk+1),');']);
%         end
%     end
%     
%     Req2=(Rtmp.*Rcu33)./(Rtmp+Rcu33)+Rsy+R2;
%     Req3=Rtemp+Rsy+R2;
%     
% tetaW1=Req1*Pjs+tetaFrame;
% tetaW2=Req2*Pjs+tetaFrame;
% tetaW3=Req3*Pjs+tetaFrame;
%%

Gcu21=1/(Rst1+Rcu1);
Gcu22=1/Rcu2;
Gcu23=1/Rcu3;
Gcu1=1/Rcu01;
Gst2=1/Rst2;
A=[Gcu21+Gcu22+Gcu23,-Gcu21-Gcu22,-Gcu23;
    -Gcu21-Gcu22,Gcu21+Gcu22+Gst2,-Gst2;
    -Gcu23,-Gst2,Gcu23+Gst2+1/(Rsy+R2)];
B=[Pjs;0;tetaFrame/(Rsy+R2)];

% AA=[Gcu1,-Gcu1,0,0;
%     -Gcu1,Gcu21+Gcu22+Gcu23+Gcu1,-Gcu21-Gcu22,-Gcu23;
%     0,-Gcu21-Gcu22,Gcu21+Gcu22+Gst2,-Gst2;
%     0,-Gcu23,-Gst2,Gcu23+Gst2+1/(Rsy+R2)];
% BB=[1/3*Pjs;2/3*Pjs;0;tetaFrame/(Rsy+R2)];
%
% tetaW44=AA\BB;

tetaW4=A\B;
Temp = tetaW4(1)+Rcu01*Pjs*Acu1/Acu;
% tetaCU2=tetaW4(1);
% Rtmp222=Rcu01;
% teta_th1=tetaW4(2);
% teta_th2=tetaW4(3);
% teta_y=teta_th2-(teta_th2-tetaFrame)*Rsy/(Rsy+R2);
% keyboard
% Temp.teta=[tetaCU;tetaW4;teta_y;tetaFrame];
% Temp.legend={'Cu/nslot layer 1','slot layer 2','node3 tooth','node4 tooth', 'yoke temperature','Frame temperature'};