% Copyright 2018
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

function [map] = syrmDesign_Vtype(dataSet)
% 
% [map] = syrmDesign_Vtype(dataSet)
% 
% Preliminary design for IPM Vtype machines. References:
% - IEMDC 2019 paper (Torino - Raleigh project)

%% inputs

if (dataSet.NumOfLayers~=1)
    error('syrmDesign_Vtype not support IPM with multiple layers!!!')
end

mu0 = 4e-7*pi;              % air permeability

Bfe = dataSet.Bfe;          % no-load steel loading (yoke flux density [T])
kj = dataSet.ThermalLoadKj; % thermal loading (copper loss per stator outer surface [W/m^2])
kt = dataSet.kt;            % kt = wt/wt_unsat
Bs = 2;                   % maximum flux density to consider the iron as ideal [T]

[~, ~, geo, per, mat] = data0(dataSet);

R          = geo.R;
p          = geo.p;
q          = geo.q;
acs        = geo.acs;
g          = geo.g;
l          = geo.l;
avv        = geo.win.avv;
kcu        = geo.win.kcu;
Ns         = geo.win.Ns;
% ns         = geo.ns;
ttd        = geo.ttd;
tta        = geo.tta;
RaccordoFC = geo.SFR;
pont0      = geo.pont0;
pontT      = geo.pontT(1);
pontR      = geo.pontR(1);
wrib       = pontT(1)+pontR(1)/2;
%alpha_pu   = geo.dalpha_pu;
hfe_min    = geo.hfe_min;
% Ar         = geo.Ar;
Nbob       = Ns/p/q/2;                                              % conductors in slot per layer
n3ph = geo.win.n3phase;

% Br = mat.LayerMag.Br;
% mur = mat.LayerMag.mu;
tempcu = per.tempcu;

if isfield(mat.LayerMag,'temp')
    tempPM = per.tempPP;
    BrTemp = mat.LayerMag.temp;
    Br = interp1(BrTemp.temp,BrTemp.Br,tempPM);
else
    Br = mat.LayerMag.Br;
end

loadpu = dataSet.CurrLoPP;                                          % current load in p.u. of i0

%flags for Vtype design

flag_beta=1;    % flag_beta=0 --> beta=90°
                % flag_beta=1 --> beta to maximize anisotropy

flag_alfa=1;    % flag_alfa=0 --> alpha_pu fixed
                % flag_alfa=1 --> alpha_pu computed from the wt                

% design domain according to x and hc/g
m = 31; n = 21;                                                     % m x n grid of evaluated machines
h = linspace(dataSet.bRange(1),dataSet.bRange(2),m);                % iron/copper split factor
x = linspace(dataSet.xRange(1),dataSet.xRange(2),n);                % rotor/stator split factor

%         [~, ~, geo, per, ~] = data0(dataSet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parametric analysis: design domain (x,hc/g)
[xx,hh] = meshgrid(x,h);

%     [kw, ~] = calcKwTh0(avv,6*p*q,p);
[kw, ~] = calcKwTh0(geo);
kw = kw(1);

r = R*xx;
rocu = 17.8*(234.5 + tempcu)/(234.5+20)*1e-9;                           % resistivity of copper [Ohm m]
ssp = r * pi/(3*p*q);                                                   % stator slot pitch (x,b)
sso = ssp * acs;                                                        % stator slot opening (x,b)
kc = ssp./(ssp-2/pi*g*(sso/g.*atan(sso/(2*g))-log(1+(sso/(2*g)).^2)));	% Carter coefficient (x,b)


cos_x0 = cos(pi/2/p);
sin_x0 = sin(pi/2/p);
kshaft = (1/cos_x0-sqrt(((1-cos_x0^2)/cos_x0)^2+sin_x0^2));
ArLim = geo.R*xx*kshaft;                                               % (max) shaft radius (x,b) [mm]
Ar = ArLim;

% Ar = Ar*ones(size(xx));
% Ar(Ar>ArLim)=ArLim(Ar>ArLim);
% Ar=ArLim;

hcMin = 2*pont0;
la    = R*xx-ArLim-hfe_min-pont0;
hcMax = la/2;

hc = hh.*g; % PMs thickness [mm]
hc_pu = hh;

%% rotor design + stator design (iron)
% check for multidimensional matrix!!!

alphaVect = linspace(0.5,0.9,101);
kalphaVect = sin(alphaVect*pi/2/p)./(alphaVect*pi/2/p);

alpha_pu = zeros(size(xx));
beta     = zeros(size(xx));
Bsq      = zeros(size(xx));
wt       = zeros(size(xx));

for rr=1:m
    for cc=1:n
        wtrVect = (r(rr,cc)-hc(rr,cc)/2-pontT)*sin((1-alphaVect)*pi/2/p)-hc(rr,cc)/2;
        if flag_beta==0
            betaVect=zeros(size(alphaVect));
        elseif flag_beta==1
            betaVect=calc_Vtype_angle_lim(geo,alphaVect*pi/2/p*180/pi,r(rr,cc)*ones(size(alphaVect)),ArLim(rr,cc)+wtrVect,hc(rr,cc)*ones(size(alphaVect)));
            betaVect=betaVect*180/pi;
        end
        
        BsqVect = Br*(kalphaVect.*(1-pontT/r(rr,cc)-hc(rr,cc)/r(rr,cc)/2)-pontR/2/r(rr,cc)./(alphaVect*pi/2/p)-sind(betaVect)./(alphaVect*pi/2/p).*(wrib/r(rr,cc)*Bs/Br+pont0/r(rr,cc)))./(sind(betaVect)+kc(rr,cc)*g/hc(rr,cc).*kalphaVect);
        
        wtVect = 2*pi*r(rr,cc)/(6*p*q*n3ph)*BsqVect/Bfe*kt;
        fObj = wtrVect-wtVect/2;
        if flag_alfa==0
            alpha_pu(rr,cc) = geo.dalpha_pu;
        else
            alpha_pu(rr,cc) = interp1(fObj,alphaVect,0);
        end
        beta(rr,cc)     = interp1(alphaVect,betaVect,alpha_pu(rr,cc));
        Bsq(rr,cc)      = interp1(alphaVect,BsqVect,alpha_pu(rr,cc));
        wt(rr,cc)       = interp1(alphaVect,wtVect,alpha_pu(rr,cc));
    end
end


alpha_pu(hc<hcMin) = NaN;
beta(hc<hcMin)     = NaN;
Bsq(hc<hcMin)      = NaN;
wt(hc<hcMin)       = NaN;

alpha_pu(hc>hcMax) = NaN;
beta(hc>hcMax)     = NaN;
Bsq(hc>hcMax)      = NaN;
wt(hc>hcMax)       = NaN;


ly = pi/2*R/p*xx.*alpha_pu.*Bsq/Bfe;
lt = R-r-g-ly;
alpha_pu(lt<2*ttd+RaccordoFC) = NaN;
beta(lt<2*ttd+RaccordoFC)     = NaN;
Bsq(lt<2*ttd+RaccordoFC)      = NaN;
wt(lt<2*ttd+RaccordoFC)       = NaN;
ly(lt<2*ttd+RaccordoFC)       = NaN;
lt(lt<2*ttd+RaccordoFC)       = NaN;

%% base inductance, d-axis parameters

kalpha = sin(alpha_pu*pi/2/p)./(alpha_pu*pi/2/p);

Lbase  = 6*mu0/pi*(kw*Ns/p)^2*(R*1e-3)*(l*1e-3)./(kc*(g*1e-3)).*xx; % base inductance (solid rotor) [H]
fM     = 2*(R/1000)*(l/1000)*kw*Ns/p*xx.*4/pi.*sin(alpha_pu*pi/2).*Bsq;
Ldpu   = 1-sin(alpha_pu*pi/2)./(1+(kc*g)./hc./sind(beta).*kalpha);



%Ldpu  = (1-sin(alpha_pu*pi/2))+kc*g./hh.*kalpha./(sind(beta)+kc*g./hc.*kalpha)-Bsq./id*pi/6/mu0*p^2/kw/Ns.*(kc*g.*sin(alpha_pu*pi/2/p))./(sind(beta)+kc*g./hc.*kalpha);

%% slot area and leakage computation + copper temperature estimation

Aslots  = zeros(m,n);
d1      = zeros(m,n);
c0      = zeros(m,n);
c1      = zeros(m,n);
c2      = zeros(m,n);
dTempCu = zeros(m,n);

for rr=1:m
    for cc=1:n
        % slot area evaluation
        geo0 = geo;
        geo0.r  = r(rr,cc);
        geo0.lt = lt(rr,cc);
        geo0.wt = wt(rr,cc);
        try
            [tmp1,tmp2] = drawSlot(geo0);
            Aslots(rr,cc) = 6*p*q*n3ph*tmp1.Aslot;
            d1(rr,cc) = tmp2.d1;
            c0(rr,cc) = tmp2.c0;
            c1(rr,cc) = tmp2.c1;
            c2(rr,cc) = tmp2.c2;
        catch
            d1(rr,cc) = NaN;
            c0(rr,cc) = NaN;
            c1(rr,cc) = NaN;
            c2(rr,cc) = NaN;
            Aslots(rr,cc) = NaN;
        end
        
        % copper overtemperature
        geo0=geo;
        geo0.r  = geo.R*xx(rr,cc);
        geo0.wt = wt(rr,cc);
        geo0.lt = lt(rr,cc);
        geo0.ly = ly(rr,cc);
        
        dTempCu(rr,cc) = temp_est_simpleMod(geo0,per);
    end
end

dTempCu=dTempCu-per.temphous;

%% slot leakage
d0 = ttd;
d2 = lt- d0 - d1;
betaSlot = c1./c2;
h = (betaSlot.^2-betaSlot.^4/4-log(betaSlot)-0.75)./((1-betaSlot).*(1-betaSlot.^2).^2);
ps = d0./c0 + d1./c0.*1./(c1./c0-1).*log(c1./c0)+d2./c2.*h;
Ls = 2*mu0*(l/1000)*Ns^2/p*ps/q/n3ph;            % According to Lipo's Book (2017)
Lspu = Ls./Lbase;

%% rated current
% Aslots = 2 * area_half_slot *6*p*q;

if q<1
    lend = 0.5*(wt+pi*(r+lt/2)*sin(pi/(6*p*q*n3ph)));
else
    lend = 2*lt+(0.5*pi*(R-ly+r)/p);                                % end turn length (x,b) [mm]
end

i0     = (kj*kcu/rocu*l./(l+lend).*pi*(R/1000).*(Aslots/1e6)/9/Ns^2).^0.5;
i0     = real(i0);
iAmp   = i0*loadpu;

J = 2*Nbob * iAmp ./ (Aslots/(q*6*p*n3ph)*kcu);                          % current density in copper [A/mm2] pk
A = 2*Nbob * iAmp ./ (r*2*pi/(q*6*p*n3ph));                              % linear current density [A/mm] pk

%% Working point + Lqpu computation
if flag_beta==0
    gamma=90;
else
    gamma=135;
end

id = iAmp*cosd(gamma);
iq = iAmp*sind(gamma);


% Ldpu  = (1-sin(alpha_pu*pi/2))+kc*g./hc.*kalpha./(sind(beta)+kc*g./hc.*kalpha);
% Ldpu  = 2*(sin(alpha_pu*pi/2).^2.*(kc*g./hc)./(sind(beta)./kalpha+kc*g./hc)+(1-sin(alpha_pu*pi/2)).^2);
% Ldpu  = alpha_pu.*sin(alpha_pu*pi/2).*(kc*g./hc)./(kc*g./hc+sind(beta)./kalpha)+(1-alpha_pu).*(1-sin(alpha_pu*pi/2));
Ldpu  = sin(alpha_pu*pi/2).*(kc*g./hc)./(kc*g./hc+sind(beta)./kalpha)+1-sin(alpha_pu*pi/2);
% Ldpu  = 8/pi^2*(sin(alpha_pu*pi/2).^2.*(kc*g./hc)./(kc*g./hc+sind(beta)./kalpha)+(1-sin(alpha_pu*pi/2)).^2./(1-alpha_pu));
ich   = fM./((Ldpu+Lspu).*Lbase);

BH=mat.Stator.BH;                   % B-H curve
muDiff=(BH(2:end,1)-BH(1:end-1,1))./(BH(2:end,2)-BH(1:end-1,2));
muDiff=[0; muDiff];
muDiff=muDiff/mu0;

Hy0   = interp1(BH(:,1),BH(:,2),Bfe);
Ht0   = interp1(BH(:,1),BH(:,2),Bfe/kt);
muy   = interp1(BH(:,1),muDiff,Bfe)*mu0;
mut   = interp1(BH(:,1),muDiff,Bfe/kt)*mu0;
Lqpu0 = 2*p*kc*g./(mu0*pi*r)./(2*p*kc*g./(mu0*pi*r)+lt./(mut.*wt)+pi./(3*p*q*muy).*(R./ly-1/2));    % Ldpu @ By==Bfe && Bt==Bfe/kt
ieq   = ich.*(Ldpu+Lspu)./(Lqpu0+Lspu);
Hy    = Hy0./ieq.*abs(iq);
Ht    = Ht0./ieq.*abs(iq);
muy   = interp1(BH(:,2),muDiff,Hy)*mu0;
mut   = interp1(BH(:,2),muDiff,Ht)*mu0;
Lqpu  = (p*kc*g)./(mu0*pi*r)./((p*kc*g)./(mu0*pi*r)+lt./(mut.*wt)+pi./(6*p*q*muy).*(R./ly-1/2));


%% Performance computation

Ld = (Ldpu+Lspu).*Lbase;
Lq = (Lqpu+Lspu).*Lbase;
fd = (fM+Ld.*id);
fq = (Lq.*iq);

T     = 3/2*p.*(fd.*iq-fq.*id);
delta = atan2(fq,fd);
PF    = sin(gamma-delta);


%% save all the results in the map structure
map.xx       = xx;
map.bb       = hh;
map.beta     = beta.*ones(size(xx))*pi/180;
map.alpha_pu = alpha_pu;
map.T        = T;
map.PF       = PF;
map.id       = id;
map.iq       = iq;
map.fd       = fd;
map.fq       = fq;
map.fM       = fM;
map.wt       = wt;
map.lt       = lt;
map.hc_pu    = hc_pu;
map.hc       = hc;
map.Ar       = Ar;
map.geo      = geo;
map.ly       = ly;
map.Lbase    = Lbase;
map.Ldpu     = Ldpu;
map.Lqpu     = Lqpu;
map.Lspu     = Lspu;
map.dTempCu  = dTempCu;
map.Aslots   = Aslots;
map.J        = J;
map.A        = A;
map.ich      = ich;
map.i0       = i0;
map.iAmp     = iAmp;
map.gamma    = atan2(map.iq,map.id)*180/pi;
map.Bsq      = Bsq;
map.flag_PM  = flag_beta;
map.dataSet  = dataSet;
                