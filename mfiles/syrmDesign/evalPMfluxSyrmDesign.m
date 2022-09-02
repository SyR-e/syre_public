% Copyright 2022
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

function [fM] = evalPMfluxSyrmDesign(tmp)
debug = 0;
nlay    = tmp.nlay;
Br      = tmp.Br;
Bs      = tmp.Bs.*ones(1,nlay);
% Bs      = [1.6 1.8 1.3]*.91;
mum     = 1;
mu0     = 4*pi*10^-7;


p       = tmp.p;
Ns      = tmp.Ns;
g       = tmp.g/1000;
sk      = tmp.sk/1000;
hc      = tmp.hc/1000;
r       = tmp.r/1000;
dalpha  = tmp.dalpha*pi/180;
l       = tmp.l/1000;
kc      = tmp.kc;
kw      = tmp.kw;
pontT   = tmp.pontT/1000;
pontR   = tmp.pontR/1000;
Tfillet = tmp.Tfillet/1000;
Rfillet = tmp.Rfillet/1000;

xPM     = sum(tmp.PMdim)./tmp.sk;
% xPM     = [1 1 1];
yPM     = ones(size(xPM));
xRibs   = (pontT+pontR/2)./sk;
% xAir    = (Tfillet+Rfillet)./sk; 
xAir    = 1-xPM; 

filt = (xAir < (Tfillet+Rfillet)./sk);
xAir(filt)    = (Tfillet(filt)+Rfillet(filt))./sk(filt);
xPM(filt) = xPM(filt) - xAir(filt);

% Magnet
% Rm      = hc./(mu0*mum*sk.*(xPM-xRibs-xAir)*l);
Rm      = hc./(mu0*mum*sk.*(xPM)*l);
Rg      = g*kc./(mu0.*dalpha*r*l);
% phim    = Br*sk.*(xPM-xRibs-xAir)*l;
phim    = Br*sk.*(xPM)*l;

% Ribs
Rr      = hc./(mu0*sk.*xRibs*l);
phir    = Bs.*sk.*xRibs*l;

% Air
% Ra      = hc./(mu0*(sk.*(1-xPM)-xRibs)*l);
Ra      = hc./(mu0*sk.*xAir*l);

% Total 
Rb      = 1./(1./Ra+1./Rm+1./Rr);
% Rb      = 1./(1./Rm+1./Rr);
phib    = phim-phir;

Rm      = fliplr(Rm);
Rg      = fliplr(Rg);
Rb      = fliplr(Rb);
phib    = fliplr(phib);

% Reluctance matrix
RR      = zeros(nlay);
FF      = zeros(1,nlay);
for ii=1:(nlay-1)
    RR(ii,ii)   =  1/Rb(ii)+1/Rg(ii)+1/Rb(ii+1);
    RR(ii,ii+1) = -1/Rb(ii+1);
    FF(ii)      =  phib(ii)-phib(ii+1);
end
for ii=2:(nlay)
    RR(ii,ii-1) = -1/Rb(ii);
end
RR(nlay,nlay)= 1/Rb(nlay)+1/Rg(nlay);
FF(nlay) = phib(nlay);


rk      = RR\FF';
phiGAP  = rk'./Rg;
phiGAP  = fliplr(phiGAP);
% fM      = 2*sum(phiGAP)*Ns*kw;

%% First harmonic factor
xmax = pi; ph = 0;
x = linspace(0,xmax,1000);
Fs = sin(x - ph);
alphaElt = tmp.alpha*p*pi/180;

angles = fliplr(alphaElt);
angles = pi/2 - angles;
angles = [0 angles];
if angles(end) < pi/2
    angles = [angles pi/2];
end
% angles = [angles(1:nlay+1) angles(end)];
levels = (cos(angles(1:end-1)) - cos(angles(2:end)))./(angles(2:end) - angles(1:end-1));
levels(1) = 0;
temp = fliplr(levels);
temp = [levels temp(2:end)];
levels = temp;
angles = angles(2:end);
angles_all = [angles(1:end-1) pi-(fliplr(angles(1:end-1)))];

y = zeros(1,length(x));
for jj = 2:length(angles_all)
    y((x < angles_all(jj)) & (x > angles_all(jj-1))) = levels(jj);
end

if (debug)
    figure()
    figSetting
    plot(x*180/pi,abs(Fs),'LineWidth',2), grid on, axis([ph 180 0 1]);
    hold on
    plot(x*180/pi,abs(y),'k','LineWidth',2)
    set(gca,'PlotBoxAspectRatio',[1 0.5 1])
    xlabel('rotor angular coordinate - elt degrees')
    ylabel('MMF - Am')
end

levels = [0 phiGAP];
levels = cumsum(levels);
temp = fliplr(levels);
temp = [levels temp(2:end)];
levels = temp/sum(phiGAP);

x = linspace(-pi,pi,1001);

xtmp = linspace(0,xmax,1000);
ytmp = zeros(1,length(xtmp));
for jj = 2:length(angles_all)
    ytmp((xtmp < angles_all(jj)) & (xtmp > angles_all(jj-1))) = levels(jj);
end

% y=zeros(1,length(x));
xtmp=xtmp-pi/2;
y=interp1(xtmp,ytmp,x);
y(x>+pi/2)=-interp1(xtmp,ytmp,fliplr(abs(x(x>+pi/2)-pi/2)));
y(x<-pi/2)=-interp1(xtmp,ytmp,fliplr(abs(x(x<-pi/2)+pi/2)));
x = x+pi/2;
if (debug)
    figure
    figSetting
    plot(x*180/pi,y,'-.g','LineWidth',1.5,'DisplayName','unsaturated staircase')
end

% saturation
y(y>+1) = +1;
y(y<-1) = -1;

% mean and first harmonic factor
kfm=mean(abs(y));
kf1=1/pi*trapz(x,y.*sin(x));

if (debug)
    plot(x*180/pi,y,'-g','LineWidth',2,'DisplayName','saturated staircase')
    plot([-180 180],+[1 1],'--k','LineWidth',1,'HandleVisibility','off')
    plot([-180 180],-[1 1],'--k','LineWidth',1,'HandleVisibility','off')
    plot(x*180/pi,kfm*ones(1,length(x)),'--r','LineWidth',1,'DisplayName','mean: y=b*kfm')
    plot(x*180/pi,kf1*sin(x),'-r','LineWidth',1,'DisplayName','1^{st} harmonics: y=b*kf1*cos(x)')
    legend('show','Location','South')
    ylim([0 1])
    xlim([0 180])
end
% kf1=1;
fM      = 2*sum(phiGAP)*Ns*kw*kf1;
