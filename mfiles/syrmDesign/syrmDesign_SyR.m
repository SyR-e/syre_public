% Copyright 2016
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

function [map] = syrmDesign_SyR(dataSet,flags)
%
% [map] = syrmDesign_SyRM(dataSet)
%
% Preliminary design for SyR machines. References:
% - Vagati's Tutorial 1994
% - ECCE 2018 paper (FEAfix)

% input
if nargin==1
    %flags for SyRM design
    flag_kw=1;      % flag_kw=0 --> use Vagati's equations, with kw=pi/(2*sqrt(3))
                    % flag_kw=1 --> use the winding factor

    flag_pb=1;      % flag_pb=0 --> hc             = constant (useful for adding PMs)
                    % flag_pb=1 --> sk/hc          = constant (reduce harmonic content)
                    % flag_pb=2 --> hc/(df*sk^0.5) = constant (reduce Lfq)

    flag_dx=1;      % flag_dx=0 --> dx=0
                    % flag_dx=1 --> constant rotor carrier width
                    % flag_dx=2 --> rotor carrier width proportional to sine integral
                    % flag_dx=3 --> rotor carrier width proportional to flux of d-axis staircase (kt needed)

    flag_ks=1;      % flag_ks=0 --> no saturation factor used
                    % flag_ks=1 --> saturation factor enabled
else
    flag_kw = flags.kw;
    flag_pb = flags.pb;
    flag_dx = flags.dx;
    flag_ks = flags.ks;
end

mu0 = 4e-7*pi;                  % air permeability
Bs = 2.0;                       % ribs flux density [T]
Bfe = dataSet.Bfe;              % steel loading (yoke flux density [T])
kj  = dataSet.ThermalLoadKj;    % thermal loading (copper loss/stator outer surface [W/m^2])
kt  = dataSet.kt;               % kt = wt/wt_unsat

[~, ~, geo, per, mat] = data0(dataSet);

R     = geo.R;
p     = geo.p;
q     = geo.q;
acs   = geo.acs;
g     = geo.g;
l     = geo.l;
avv   = geo.win.avv;
kcu   = geo.win.kcu;
Ns    = geo.win.Ns;
% ns    = geo.ns;
ttd   = geo.ttd;
tta   = geo.tta;
RaccordoFC = geo.SFR;
nlay  = geo.nlay;
pont0 = geo.pont0;
pont  = geo.pont0;%+min(geo.pont);
n3ph  = geo.win.n3phase;

tempcu = per.tempcu;

% design domain according to b and x
m = 31; n = 21;                                                     % m x n grid of evaluated machines
b = linspace(dataSet.bRange(1),dataSet.bRange(2),m);                % iron/copper split factor
x = linspace(dataSet.xRange(1),dataSet.xRange(2),n);                % rotor/stator split factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parametric analysis: design domain (x,b)
[xx,bb] = meshgrid(x,b);

[xGap,yGap,kf1,kfm] = evalBgapSyrmDesign(q,kt);                     % airgap induction shape, first harmonic factor, mean value factor

r = R*xx;                                                           % rotor radius [m]
rocu = 17.8*(234.5 + tempcu)/(234.5+20)*1e-9;                       % resistivity of copper [Ohm m]
ssp = r * pi/(3*p*q);                                               % stator slot pitch (x,b)
sso = ssp * acs;                                                    % stator slot opening (x,b)
kc = ssp./(ssp-2/pi*g*(sso/g.*atan(sso/(2*g))-log(1+(sso/(2*g)).^2)));   % Carter coefficient (x,b)
ly = pi/2*R/p*xx.*bb*kfm;                                           % yoke or back iron (x,b) [mm]
%ly = R/p*xx.*bb;                                                   % yoke [mm], do not depend on kt
wt = 2*pi*R/(6*p*q*n3ph)*xx.*bb.*kt;                                     % tooth width (x,b) [mm]
cos_x0 = cos(pi/2/p);
sin_x0 = sin(pi/2/p);
Ar = geo.R*xx*(1/cos_x0-sqrt(((1-cos_x0^2)/cos_x0)^2+sin_x0^2));    % (max) shaft radius (x,b) [mm]

lt = R*(1-xx)-g-ly;
lt(lt<2*ttd) = NaN;

% total insulation la. See readme\2015 11 03 - define la in syreDesign.pptx
% 0) la(x,b) = r - Ar - ly (total radial space r - Ar: too much insulation)
% la0 = R * xx - Ar - ly;
% 1) la(x,b) = x0 - rbeta - ly (reduced radial space x0-rbeta-Ar: too much iron)
% la1 = x0 - rbeta - Ar - ly;
x0 = r /cos(pi/2/p);                                                % center of barriers circles
geo.dalpha = geo.dalpha_pu*(90/p);                                  % [mec degrees]
beta_temp = atand(r*sind(geo.dalpha(1))./(x0 - r * cosd(geo.dalpha(1))));
rbeta = (x0 - r*cosd(geo.dalpha(1)))./(cosd(beta_temp));            % radius of barrier circle
% 2) 1st carrier takes 1-cos(p*alpha1) p.u. flux
la = (x0 - rbeta - Ar - ly*cosd(p*geo.dalpha(1)))*nlay/(nlay-0.5);
%%
% rotor design + slot evaluation + Lfqpu + Lcqpu
alpha = cumsum(geo.dalpha);                                         % alpha in syre coordinates (mech deg, zero is axis Q)
[df,da] = staircaseAnyAlpha(alpha*p*pi/180);                        % rotor staircase: set of rotor slots positions
f = cumsum((df));                                                   % stator MMF staircase
sumDf2r = sum(df.^2);
Lcqpu = 1-4/pi*sum(f.^2.*da);

alpha = cumsum(da);                                                 % alpha defined as in Vagati tutorial (0 = d axis)
geo.alpha = cumsum(geo.dalpha);                                     % alpha in syre coordinates (mech deg, zero is axis Q)

% initializing matrix
% q-axis
Lfqpu = zeros(m,n);
% rotor
hc = cell(m,n);
sk = cell(m,n);
hf = cell(m,n);
dx = cell(m,n);
hc_pu = cell(m,n);
Bx0 = cell(m,n);
delta = zeros(m,n);
lr = zeros(m,n);
% slot
Aslots = zeros(m,n);
d1 = zeros(m,n);
c0 = zeros(m,n);
c1 = zeros(m,n);
c2 = zeros(m,n);
% copper temperature
dTempCu = zeros(m,n);
% additional ribs
pontR = zeros(m,n);

dfQ=fliplr(df); % df using syre conventions (alpha=0 is the q-axis)

% design of the single machines to get slot area, barriers dimensions Lfqpu and copper overtemperature
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
        
        % rotor design
        % sk
        beta = calc_apertura_cerchio(geo.alpha*pi/180,geo.R*xx(rr,cc),x0(rr,cc));
        rbeta = (x0(rr,cc) - R*xx(rr,cc) * cosd(geo.alpha))./(cos(beta));
        [xpont,ypont] = calc_intersezione_cerchi(R*xx(rr,cc)-geo.pont0, rbeta, x0(rr,cc));
        if strcmp(dataSet.TypeOfRotor,'Circular')
            sk{rr,cc}=rbeta.*beta;
            lr(rr,cc)=mean(sk{rr,cc}(end-1:end)); % length of the inner flux carrier (for saturation factor)
        elseif strcmp(dataSet.TypeOfRotor,'Seg')|| strcmp(dataSet.TypeOfRotor,'ISeg')
            rpont_x0=sqrt(ypont.^2+(x0(rr,cc)-xpont).^2);
            [alphapont,rpont] = cart2pol(xpont,ypont);
            Bx0=x0(rr,cc)-(rpont_x0);
            mo=1;
            y=ypont+mo*(Bx0-xpont);
            xBmk=Bx0;
            yBmk=y;
            skv{rr,cc}=calc_distanza_punti_altern(xBmk,yBmk,Bx0,zeros(size(Bx0)));
            sko{rr,cc}=calc_distanza_punti_altern(xpont,ypont,xBmk,yBmk);
            lrv(rr,cc)=mean(skv{rr,cc}(end-1:end));
            lro(rr,cc)=mean(skv{rr,cc}(end-1:end));

            sk{rr,cc}=skv{rr,cc}+sko{rr,cc};
            lr(rr,cc)=lrv(rr,cc)+lro(rr,cc);
        end
        
        clear Bx0
        switch flag_pb
            case 0 % hc = cost
                hc{rr,cc}=la(rr,cc)/geo.nlay.*ones(1,geo.nlay);
                %disp('flux barrier design: hc = cost')
            case 1 % pb = cost
                hc{rr,cc}=la(rr,cc)/sum(sk{rr,cc}).*sk{rr,cc};
                %disp('flux barrier design: pbk = hc/sk = cost')
            case 2 % min Lfq
                hc{rr,cc}=la(rr,cc)/sum(dfQ.*sk{rr,cc}.^0.5).*(dfQ.*sk{rr,cc}.^0.5);
                %disp('flux barrier desig: hc/(df*sk) = cost')
        end
        
        rpont_x0=sqrt(ypont.^2+(x0(rr,cc)-xpont).^2);
        Bx0{rr,cc}=x0(rr,cc)-(rpont_x0);
        hc_min=(R*xx(rr,cc)-Ar(rr,cc)-(geo.R-geo.R*xx(rr,cc)-geo.g-lt(rr,cc)))/geo.nlay/4;
        hfe_min=2*geo.pont0;
        delta(rr,cc)=(0.5*hc{rr,cc}(1)+sum(hc{rr,cc}(2:end-1))+0.5*hc{rr,cc}(end)-hc_min*(geo.nlay-1))/(Bx0{rr,cc}(1)-Bx0{rr,cc}(end)-hfe_min*(geo.nlay-1)-hc_min*(geo.nlay-1));
        hc_pu{rr,cc}=hc{rr,cc}*(delta(rr,cc)*geo.nlay)/sum(hc{rr,cc});
        
        % dx (flux carrier design)
        alphad=[0 90-fliplr(geo.alpha)*geo.p 90];                   % 0<=alphad<=90 [° elt]
        r_all=geo.R*xx(rr,cc);
        for ii=1:geo.nlay
            r_all=[r_all Bx0{rr,cc}(ii)+hc{rr,cc}(ii)/2 Bx0{rr,cc}(ii)-hc{rr,cc}(ii)/2];
        end
        r_all=[r_all Ar(rr,cc)];
        hf0=r_all(1:2:end)-r_all(2:2:end);
        
        switch flag_dx
            case 0
                hf{rr,cc}=hf0;
            case 1 % constant iron
                hf_cost=[ones(1,geo.nlay-1) 0.5]/(geo.nlay-0.5);
                hf{rr,cc}=hf_cost*sum(hf0(2:end));
                %disp('flux carrier design: Fe = cost')
            case 2 % iron proportional to first harmonic flux
                level=(cosd(alphad(1:end-1))+cosd(alphad(2:end)))/2.*(alphad(2:end)-alphad(1:end-1));
                level_pu=level/sum(level);
                hfTemp=fliplr(level_pu)*sum(hf0);
                hf{rr,cc}=hfTemp(2:end);
                %disp('flux carrier design: Fe proportional to first harmonic flux')
            case 3 % iron proportional to flux
                level_pu=evalSatStairCase(xGap,yGap,alphad);
                hf{rr,cc}=fliplr(level_pu)*sum(hf0);
                %disp('flux carrier design: Fe proportional to flux')
        end
        
        for ii=geo.nlay:-1:1
            if ii==geo.nlay
                B1tmp=Ar(rr,cc)+hf{rr,cc}(ii);
            else
                B1tmp=B2k(ii+1)+hf{rr,cc}(ii);
            end
            dxtmp=1-2/hc{rr,cc}(ii)*(Bx0{rr,cc}(ii)-B1tmp);
            B1ktmp=Bx0{rr,cc}(ii)-hc{rr,cc}(ii)/2+dxtmp*hc{rr,cc}(ii)/2;
            B2ktmp=Bx0{rr,cc}(ii)+hc{rr,cc}(ii)/2+dxtmp*hc{rr,cc}(ii)/2;
            if B1ktmp>(Bx0{rr,cc}(ii)-geo.pont0)
                B1new=Bx0{rr,cc}(ii)-geo.pont0;
                dxtmp=1-2/hc{rr,cc}(ii)*(Bx0{rr,cc}(ii)-B1new);
            elseif B2ktmp<(Bx0{rr,cc}(ii)+geo.pont0)
                B2new=Bx0{rr,cc}(ii)+geo.pont0;
                dxtmp=2/hc{rr,cc}(ii)*(B2new-Bx0{rr,cc}(ii))-1;
            end
            dx{rr,cc}(ii)=dxtmp;
            B1k(ii)=Bx0{rr,cc}(ii)-hc{rr,cc}(ii)/2+dx{rr,cc}(ii)*hc{rr,cc}(ii)/2;
            B2k(ii)=Bx0{rr,cc}(ii)+hc{rr,cc}(ii)/2+dx{rr,cc}(ii)*hc{rr,cc}(ii)/2;
        end
        
        if flag_dx==0
            dx{rr,cc}=zeros(1,geo.nlay);
            %disp('flux carrier design: dx=0')
            hf{rr,cc}=hf0;
        end
        
        % q-axis flow-through inductance [pu]
        Lfqpu(rr,cc) = 4/pi*p*g*kc(rr,cc)/(xx(rr,cc)*R)*(sum((dfQ).^2.*(sk{rr,cc}./hc{rr,cc})));
        
        % copper overtemperature
        geo0=geo;
        geo0.r  = geo.R*xx(rr,cc);
        geo0.wt = wt(rr,cc);
        geo0.lt = lt(rr,cc);
        geo0.ly = ly(rr,cc);
        
        dTempCu(rr,cc) = temp_est_simpleMod(geo0,per);
        
        % radial ribs (high speed motors)
        temp.B1k=B1k;
        temp.B2k=B2k;
        temp.Bx0=Bx0;
        
        geo0.x0 = x0(rr,cc);
        geo0.r = r(rr,cc);
        geo0.RotType = 'Circular';
        geo0.hc = hc{rr,cc};
        if rr==1&&cc==1
            warning('Ribs computed for circular geometry');
        end
        
        [~,geo0] = calc_ribs_rad_fun(geo0,mat,temp);
        
        pontR(rr,cc) = mean(geo0.pontR);
        
    end
end

dTempCu=dTempCu-per.temphous;

%%
% d axis
if flag_kw
    [kw, ~] = calcKwTh0(geo);
    kw = kw(1);
else
    kw = pi/2/sqrt(3);
    warning('Winding factor not evaluated')
end

Ht = interp1(mat.Stator.BH(:,1),mat.Stator.BH(:,2),Bfe./kt);
Hy = interp1(mat.Stator.BH(:,1),mat.Stator.BH(:,2),Bfe);
if flag_ks
    %ksat = 1+mu0*pi/2*(Hy*pi/(6*p*q)*(R-ly/2) + Ht*lt + Hy*lr)./(bb*Bfe*kf1.*kc*g);
    ksat = 1+mu0*(Hy*pi/(6*p*q*n3ph)*(R-ly/2) + Ht*lt + Hy*lr)./(bb*Bfe*kf1.*kc*g);
else
    ksat=ones(size(xx));
    warning('saturation factor ksat not evaluated!!!')
end

Fmd = 2*(R*1e-3)*(l*1e-3)*kw*Ns*Bfe.*kf1/p.*xx.*bb;                 % flux linkage [Vs]
id = pi*Bfe*kc*(g*1e-3)*p.*ksat/(mu0*3*kw*Ns).*bb;                    % id [A]
%id = 2*Bfe*kc*(g*1e-3)*p.*ks/(mu0*3*kw*Ns).*bb;                    % id [A]

Lmd = Fmd./id;                                                      % d magnetization inductance Lmd [Vs]

Lbase=Lmd.*ksat;% base inductance for analytical model (unsaturated)
%Lbase=mu0*pi*(R*1e-3)*(l*1e-3)*Ns^2./(2*p^2.*kc*(g*1e-3)).*x;

% stator design
if q<1
    lend = 0.5*(wt+pi*(r+lt/2)*sin(pi/(6*p*q*n3ph)));
else
    lend = 2*lt+(0.5*pi*(R-ly+r)/p);                                % end turn length (x,b) [mm]
end

% calculate slot area and rated current
% Aslots = 2 * area_half_slot *6*p*q;

%kj = Loss/(2*pi*R*l)*1e6;                                           % specific loss (x,b) [W/m2]
K = sqrt(kcu*kj/rocu*l./(l+lend));                                  % factor K [] (x,b)

i0 = pi/(3*Ns)*(R/1000)^1.5*K.*sqrt(Aslots/(pi*R^2));               % rated current i0 [A] pk
i0=real(i0);
loadpu = dataSet.CurrLoPP;                                          % current load in p.u. of i0

id(id>loadpu*i0)=loadpu*i0(id>loadpu*i0);
gamma=acos(id./(loadpu*i0));
iq=loadpu*i0.*sin(gamma);                                           % q-axis current [A] pk

%iq = sqrt((loadpu*i0).^2 - id.^2); iq = real(iq);                  % q-axis current [A] pk
Am = Aslots.* id./(loadpu*i0);                                      % slots area dedicated to id [mm2]

Nbob  = Ns/p/q/2;                                                   % conductors in slot per layer
J = 2*Nbob * loadpu*i0 ./ (Aslots/(q*6*p)*kcu);                     % current density in copper [A/mm2] pk
A = 2*Nbob * loadpu*i0 ./ (r*2*pi/(q*6*p));                         % linear current density [A/mm] pk

% tangential ribs effect
Lrib = 4/pi*kw*Ns*((2*pont+pontR)*1e-3)*(l*1e-3)*Bs./iq;

% q axis
Lmq=Lbase.*(Lcqpu+Lfqpu)+Lrib;

kdq=1-Lmq./Lmd;

T = 3/2*p*(kdq.*Fmd.*iq);

% stator leakage inductance Ls
[dfs] = staircaseRegular(6*q);                                      % stator staircase
f = cumsum(dfs);
sumDf2s = sum(dfs.^2);
d0 = ttd;
d2 = lt- d0 - d1;
beta = c1./c2;
h = (beta.^2-beta.^4/4-log(beta)-0.75)./((1-beta).*(1-beta.^2).^2);
ps = d0./c0 + d1./c0.*1./(c1./c0-1).*log(c1./c0)+d2./c2.*h;

% Lspu = 4/pi*p*kc*g*sumDf2s.*ps./r;                                  % Ls/Lmd: p.u. leakage inductance

% Ls = 2*mu0*(l/1000)*Ns^2/p*ps.*sumDf2s;     % According to Vagati's Tutorial (1994)
Ls = 2*mu0*(l/1000)*Ns^2/p*ps/q;            % According to Lipo's Book (2017)
Lspu = Ls./Lbase;

% power factor
Ld = Lmd+Ls;
Lq = Lmq+Ls;

csi = Ld./Lq;
gamma = atand(iq./id);                                              % current phase angle [deg]
delta = atand((Lq.*iq)./(Ld.*id));                                  % flux linkage phase angle [deg]

PF = sind(gamma-delta);                                             % PF @ gamma (same gamma as torque)
PFmax = (Ld-Lq)./(Ld+Lq);                                           % PF @ max PF gamma
% PF2 = (csi-1)./sqrt(csi.^2.*(sind(gamma)).^-2+(cosd(gamma)).^-2);  % alternative formula

%% save all the results in the map structure
map.xx      = xx;
map.bb      = bb;
map.T       = T;
map.PF      = PF;
map.id      = id;
map.iq      = iq;
map.fd      = Ld.*id;
map.fq      = Lq.*iq;
map.wt      = wt;
map.lt      = lt;
map.la      = la;
map.Ar      = Ar;
map.hc_pu   = hc_pu;
map.dx      = dx;
map.geo     = geo;
map.ly      = ly;
map.ksat    = ksat;
map.Lbase   = Lbase;
map.Lmd     = Lmd;
map.Lcqpu   = Lcqpu;
map.Lfqpu   = Lfqpu;
map.Lspu    = Lspu;
map.Lrpu    = Lrib./Lbase;
map.kdq     = kdq;
map.dTempCu = dTempCu;
map.Aslots  = Aslots;
map.flag_pb = flag_pb;
map.flag_dx = flag_dx;
map.A       = A;
map.J       = J;
map.kf1     = kf1*ones(size(xx));
map.kfm     = kfm*ones(size(xx));
map.iAmp    = i0*loadpu;
map.i0      = i0;
map.gamma   = atan2(map.iq,map.id)*180/pi;
map.dataSet = dataSet;
