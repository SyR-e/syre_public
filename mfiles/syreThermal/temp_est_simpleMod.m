% Soluzione rete semplificata, riferimento formule articolo Simplified
% [Boglietti et all,Thermal model for variable speed self cooled industrial
% Induction Motor]
%

function [Temp]=temp_est_simpleMod(geo,per,Tframe)
% Thermal material initialization
Thermal_param;
% Data Initialization
Pjs=per.Loss;
% Pjs=500;
tetaFrame=per.temphous;
% max speed or rated speed
omega=geo.nmax;     %rpm
% geometrical data
p=geo.p;
q=geo.q;
Kcu=geo.win.kcu;
g=geo.g*10^-3;
R=geo.R*10^-3;
L=geo.l*10^-3;
wt=geo.wt*10^-3;
lt=geo.lt*10^-3;
ly = (geo.R - geo.r - geo.lt)*1e-3;
x=geo.r/geo.R;
% ro=0.017241*10^-6;
% Bfe=1.5;
ps=2*pi*x*R/(6*p*q*geo.win.n3phase);
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
%
lN=0.0002;
ws=ps-wt;
Rse=x*R+lt;
% Aslot=(pi*R^2*(1-b*x/p).*(1+x-2*b*kt*x))/(6*p*q)
% Aslot=pi*lt.*(2*R-lt-ly)./(6*p*q*geo.win.n3phase)-wt.*lt;
if isfield(geo,'Aslot')
    Aslot=geo.Aslot*1e-6;%/(6*pi*q);
else
    Aslot=pi*lt.*(2*R-lt-ly)./(3*p*q*geo.win.n3phase)-wt.*lt;
end

Acu=Aslot*Kcu;

Rsa=x*R+g;
Rsb=Rsa+lt;
wsEq=0.5*(Rsa+Rsb)*pi/3/p/q-wt;
hsEq=Aslot./wsEq;
% hsEq=lt;
% wsEq=Aslot/lt;

wsp=wsEq.*(1-Kcu);

Ps=2*hsEq+wsEq;

teq=((1-Kcu)*Aslot./Ps);
Rcu_ir=teq./lamdaEQ./(Ps*L)/(6*p*q*geo.win.n3phase);
%
n_element=2;
for kk=1:n_element
    eval(['Rcu',num2str(kk),'=0.5*teq./lamdaEQ./(2*hsEq/',num2str(n_element),'*L)/(6*p*q*geo.win.n3phase);']);
end
eval(['Rcu',num2str(n_element+1),'=0.5*teq./lamdaEQ./(wsEq*L)/(6*p*q*geo.win.n3phase);']);
%
% Rcu11=teq./lamdaEQ./(hsEq/2*L)/(6*p*q);
% Rcu22=Rcu11;
% Rcu33=teq./lamdaEQ./(wsEq*L)/(6*p*q);
%
w=(wsEq-2*teq)/3;
h1=hsEq-teq-w;
Acu1=w*h1;
Ps1=2*hsEq-w;
Rcu01=0.5*teq/(lamdaEQ*Ps1*L)/(6*p*q*geo.win.n3phase);
% yoke
Rsy=1/(2*pi*LambdaFE*L)*log(R./Rsb);
% teeth
pir=(6*p*q)*wt.*lt./((6*p*q)*wt.*lt+(6*p*q*geo.win.n3phase)*Aslot);
Rst=1./(2*pi*LambdaFE*L*pir).*log(Rsb./Rsa);

rm=(Rsa+Rsb)/2;
fp=2*pi/(6*p*q*geo.win.n3phase);
fe=wt./Rsa;

for kk=1:n_element
    eval(['Rst',num2str(kk),'=1./(2*pi*LambdaFE*L*pir).*log((Rsa+lt*',num2str(kk/n_element),')./(Rsa+lt*',num2str((kk-1)/n_element),'));']);
end

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
[Nut,NuEw] = NusseltNumb(0);
hag = Nut*lambdaAIR/g;
Rsag=1./(Aist.*hag);
% Rrag
Rrag=1./(2*pi*x.*R*L*hag);
% Rshaft
% Rshaft=0.0645;
% rotor +shaft
Rrshf=1./(2*pi*LambdaFE*L)*log(1.5*x*R./(x*R-ly))+0.25*0.5*L./(LambdaFE*pi*0.08^2)+0.5*0.5*1.3*L./(LambdaFE*pi*0.08^2);

%% %%%%%
%% Solved network linear system
%% %%%%%
Gcu21=1/(Rst1+Rcu1);
Gcu22=1/Rcu2;
Gcu23=1/Rcu3;
Gcu1=1/Rcu01;
Gst2=1/Rst2;
A=[Gcu21+Gcu22+Gcu23,-Gcu21-Gcu22,-Gcu23;
    -Gcu21-Gcu22,Gcu21+Gcu22+Gst2,-Gst2;
    -Gcu23,-Gst2,Gcu23+Gst2+1/(Rsy+R2)];
B=[Pjs;0;tetaFrame/(Rsy+R2)];

tetaW4=A\B;
tetaCU=tetaW4(1)+Rcu01*Pjs*Acu1/Acu;
tetaCU2=tetaW4(1);
Rtmp222=Rcu01;
teta_th1=tetaW4(2);
teta_th2=tetaW4(3);
teta_y=teta_th2-(teta_th2-tetaFrame)*Rsy/(Rsy+R2);

Temp.teta=[tetaCU;tetaW4;teta_y;tetaFrame];
Temp.legend={'Cu/nslot layer 1','slot layer 2','node3 tooth','node4 tooth', 'yoke temperature','Frame temperature'};


Temp = tetaW4(1)+Rcu01*Pjs*Acu1/Acu;
