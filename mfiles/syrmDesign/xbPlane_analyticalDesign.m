% Copyright 2023
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

function [map] = xbPlane_analyticalDesign(dataSet,flags)
%
% [map] = xbPlane_analyticalDesign(dataSet,flags), previous syrmDesign_SyR and syrmDesign_SPM
%
% Preliminary design for SyR, PM-SyR, IPM and SPM machines. References:
% - Vagati's Tutorial 1994
% - ECCE 2017 paper (SPM old)
% - ECCE 2018 paper (FEAfix SyR)
% - TIA 2020 paper (FEAfix SyR)
% - ICEM 2022 paper (IPM)
% - IEMDC 2023 paper (IPM, PM-SyR)

% input
if nargin==1
    %flags for SyRM design
    f = dataSet.syrmDesignFlag;
    flag_kw=1;      % flag_kw=0 --> use Vagati's equations, with kw=pi/(2*sqrt(3))
                    % flag_kw=1 --> use the winding factor

    flag_pb=f.hc;   % flag_pb=0 --> hc             = constant (useful for adding PMs)
                    % flag_pb=1 --> sk/hc          = constant (reduce harmonic content)
                    % flag_pb=2 --> hc/(df*sk^0.5) = constant (reduce Lfq)

    flag_dx=f.dx;   % flag_dx=0 --> dx=0
                    % flag_dx=1 --> constant rotor carrier width
                    % flag_dx=2 --> rotor carrier width proportional to sine integral
                    % flag_dx=3 --> rotor carrier width proportional to flux of d-axis staircase (kt needed)

    flag_ks=f.ks;   % flag_ks=0 --> no saturation factor used
                    % flag_ks=1 --> saturation factor enabled
    
    flag_i0=f.i0;   % flag_i0=0 --> kj=constant
                    % flag_i0=1 --> J=constant
                    % flag_i0=2 --> I=constant
else
    flag_kw = flags.kw;
    flag_pb = flags.pb;
    flag_dx = flags.dx;
    flag_ks = flags.ks;
end

gammaMin = 5*pi/180;                 % minimum current angle [elt deg] to avoid iq=0

mu0      = 4e-7*pi;                  % air permeability
Bs       = 2.0;                      % ribs flux density [T]
Bfe      = dataSet.Bfe;              % steel loading (yoke flux density [T])
kj       = dataSet.ThermalLoadKj;    % thermal loading (copper loss/stator outer surface [W/m^2])
kt       = dataSet.kt;               % kt = wt/wt_unsat
J        = dataSet.CurrentDensity;   % J= slot current density, expressed in [Apk/mm^2]
% i0       = dataSet.CurrLoPP;
loadpu   = 1;                        % current load in p.u. of i0
kyr      = dataSet.RotorYokeFactor;  % kyr = lyr/lys = increase factor for rotor yoke length compared to stator yoke
ky       = dataSet.StatorYokeFactor; % ky = ly/ly_ideal
PMdimPU  = dataSet.PMdimPU;          % per-unit PM dimension
kPM      = dataSet.kPM;              % PM filling factor
Br       = dataSet.Br;               % PM remanence @ design temperature
tempPM   = dataSet.PMtemp;           % PM temperature
alpha_pu = dataSet.ALPHApu;          % alpha_pu

if Br~=0
    mat = material_properties_layer(dataSet.FluxBarrierMaterial);
    Br  = interp1(mat.temp.temp,mat.temp.Br,tempPM);
    Bd  = interp1(mat.temp.temp,mat.temp.Bd,tempPM);
    dataSet.Br = Br;
end


PMdimPU = PMdimPU./PMdimPU;
PMdimPU(isnan(PMdimPU)) = 0;
PMdimPU = kPM*PMdimPU;

[~, ~, geo, per, mat] = data0(dataSet);
% i0 = per.i0;

if flag_kw
    [kw, ~] = calcKwTh0(geo);
    kw = kw(1);
else
    kw = pi/2/sqrt(3);
    warning('Winding factor not evaluated')
end


R          = geo.R;
p          = geo.p;
q          = geo.q;
acs        = geo.acs;
g          = geo.g;
l          = geo.l;
avv        = geo.win.avv;
kcu        = geo.win.kcu;
Ns         = geo.win.Ns;
ttd        = geo.ttd;
tta        = geo.tta;
RaccordoFC = geo.SFR;
nlay       = geo.nlay;
pont0      = geo.pont0;
pont       = geo.pont0;%+min(geo.pont);
n3ph       = geo.win.n3phase;
pontT      = geo.pontT;
Nbob       = Ns/p/q/2;                                                   % conductors in slot per layer

tempcu = per.tempcu;

if length(dataSet.xRange)>1
    % design domain according to b and x
    m = 31; n = 21;                                                     % m x n grid of evaluated machines
    b = linspace(dataSet.bRange(1),dataSet.bRange(2),m);                % iron/copper split factor
    x = linspace(dataSet.xRange(1),dataSet.xRange(2),n);                % rotor/stator split factor
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % parametric analysis: design domain (x,b)
    [xx,bb] = meshgrid(x,b);
else
    xx = dataSet.xRange;
    bb = dataSet.bRange;
    m  = 1;
    n  = 1;
end

% if q>=1
%     [~,~,kf1,kfm] = evalBgapSyrmDesign(q,kt);                     % airgap flux density shape, first harmonic factor, mean value factor
% else
%     kf1 = 2/sqrt(3);    % ratio between peak arigap flux density and peak fundamental
%     kfm = 1;
% end


r = R*xx;                                                               % rotor radius [m]
rocu = 17.8*(234.5 + tempcu)/(234.5+20)*1e-9;                           % resistivity of copper [Ohm m]
ssp = r * pi/(3*p*q*n3ph);                                              % stator slot pitch (x,b)
sso = ssp * acs;                                                        % stator slot opening (x,b)
kc = ssp./(ssp-2/pi*g*(sso/g.*atan(sso/(2*g))-log(1+(sso/(2*g)).^2)));  % Carter coefficient (x,b)
% ly = pi/2*R/p*xx.*bb*kfm;                                              % yoke or back iron (x,b) [mm]
%ly = R/p*xx.*bb;                                                       % yoke [mm], do not depend on kt
% ly = 4/9*(1+sqrt(3))*R/p*xx.*bb*ky;                                     % yoke or back iron (x,b) [mm]
% ly = (sqrt(3)-pi/6)*R/p*xx.*bb*ky;                                     % yoke or back iron (x,b) [mm], elmo
ly = R/p*xx.*bb*ky;                                                       % yoke [mm], do not depend on kt, sin
wt = 2*pi*R/(6*p*q*n3ph)*xx.*bb.*kt;                                    % tooth width (x,b) [mm]
cos_x0 = cos(pi/2/p);
sin_x0 = sin(pi/2/p);
Ar = geo.R*xx*(1/cos_x0-sqrt(((1-cos_x0^2)/cos_x0)^2+sin_x0^2));    % (max) shaft radius (x,b) [mm]
ArLim = Ar;
Ar(Ar>geo.Ar) = geo.Ar;

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
la = (x0 - rbeta - ArLim - ly*kyr*cosd(p*geo.dalpha(1)))*nlay/(nlay-0.5);
%%
% rotor design + slot evaluation + Lfqpu + Lcqpu
alpha = cumsum(geo.dalpha);                                         % alpha in syre coordinates (mech deg, zero is axis Q)
[df,da] = staircaseAnyAlpha(alpha*p*pi/180);                        % rotor staircase: set of rotor slots positions
f = cumsum((df));                                                   % stator MMF staircase
sumDf2r = sum(df.^2);
Lcqpu = 1-4/pi*sum(f.^2.*da);

% alpha = cumsum(da);                                                 % alpha defined as in Vagati tutorial (0 = d axis)
geo.alpha = cumsum(geo.dalpha);                                     % alpha in syre coordinates (mech deg, zero is axis Q)

% initializing matrix
% q-axis
Lfqpu = zeros(m,n);
% rotor
hc    = cell(m,n);
sk    = cell(m,n);
hf    = cell(m,n);
dx    = cell(m,n);
hc_pu = cell(m,n);
Bx0   = cell(m,n);
skv   = cell(m,n);
sko   = cell(m,n);
lrv   = zeros(m,n);
lro   = zeros(m,n);
delta = zeros(m,n);
lr    = zeros(m,n);
% slot
Aslots = zeros(m,n);
d1 = zeros(m,n);
c0 = zeros(m,n);
c1 = zeros(m,n);
c2 = zeros(m,n);
ws = zeros(m,n);
% additional ribs
pontR   = zeros(m,n);
BrPrime = cell(m,n);
% PM remanence and flux linkage
BrAvg = zeros(m,n);
BrMax = zeros(m,n);
BrMin = zeros(m,n);
fM    = zeros(m,n);
PMdim = cell(m,n);

hcMin = zeros(m,n); % minimum barrier thickness [mm]
skMax = zeros(m,n); % inner (half) barrier length [mm]

Abar = zeros(m,n);  %total barriers area (estimation for PMs volume)

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
            ws(rr,cc) = tmp1.st;
        catch
            d1(rr,cc) = NaN;
            c0(rr,cc) = NaN;
            c1(rr,cc) = NaN;
            c2(rr,cc) = NaN;
            Aslots(rr,cc) = NaN;
            ws(rr,cc) = NaN;
        end
        
        % rotor design
        % sk
        beta = calc_apertura_cerchio(geo.alpha*pi/180,geo.R*xx(rr,cc),x0(rr,cc));
        rbeta = (x0(rr,cc) - R*xx(rr,cc) * cosd(geo.alpha))./(cos(beta));
        [xpont,ypont] = calc_intersezione_cerchi(R*xx(rr,cc)-geo.pont0, rbeta, x0(rr,cc));
        if strcmp(dataSet.TypeOfRotor,'Circular')
            sk{rr,cc}    = rbeta.*beta;
            lr(rr,cc)    = mean(sk{rr,cc}(end-1:end)); % length of the inner flux carrier (for saturation factor)
            PMdim{rr,cc} = [sk{rr,cc};zeros(size(sk{rr,cc}))].*PMdimPU;
        elseif strcmp(dataSet.TypeOfRotor,'Seg')|| strcmp(dataSet.TypeOfRotor,'ISeg')
            rpont_x0=sqrt(ypont.^2+(x0(rr,cc)-xpont).^2);
            Bx0=x0(rr,cc)-(rpont_x0);
            mo=1;
            y=ypont+mo*(Bx0-xpont);
            y = y-geo.hcShrink.*y;
            xBmk=Bx0;
            yBmk=y;
            skv{rr,cc}=calc_distanza_punti_altern(xBmk,yBmk,Bx0,zeros(size(Bx0)));
            sko{rr,cc}=calc_distanza_punti_altern(xpont,ypont,xBmk,yBmk);
            if length(skv{rr,cc})>1
                lrv(rr,cc)=mean(skv{rr,cc}(end-1:end));
            else
                lrv(rr,cc)=skv{rr,cc};
            end
            if length(sko{rr,cc})>1
                lro(rr,cc)=mean(sko{rr,cc}(end-1:end));
            else
                lro(rr,cc)=sko{rr,cc};
            end

            sk{rr,cc}    = skv{rr,cc}+sko{rr,cc};
            lr(rr,cc)    = lrv(rr,cc)+lro(rr,cc);
            PMdim{rr,cc} = [skv{rr,cc};sko{rr,cc}].*PMdimPU;
        else
            sk{rr,cc} = 0;
            PMdim{rr,cc} = 0;
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

        hcMin(rr,cc) = min(hc{rr,cc});
        skMax(rr,cc) = sk{rr,cc}(end);

        rpont_x0=sqrt(ypont.^2+(x0(rr,cc)-xpont).^2);
        Bx0{rr,cc}=x0(rr,cc)-(rpont_x0);
        hc_min=(R*xx(rr,cc)-ArLim(rr,cc)-(geo.R-geo.R*xx(rr,cc)-geo.g-lt(rr,cc)))/geo.nlay/4;
        hfe_min=2*geo.pont0;
        
        if (nlay==1)
            hc_half_min = la(rr,cc)/nlay/8;
            % max hc according to alpha min
            hc_half_max1 = (alpha*pi/180/(1+alpha*pi/180)*(r(rr,cc)-pont0));
            % (needs division by 2 .. don't know why but it works)
            hc_half_max1 = hc_half_max1 * 2;

            % max hc according to alpha max (27 Jan 2011)
            temp_alpha_hfemin = hfe_min/r(rr,cc); % rad
            temp_alpha_hc_2 = pi/(2*p) - alpha*pi/180 - temp_alpha_hfemin;
            hc_half_max2 = (temp_alpha_hc_2/(1+temp_alpha_hc_2)*(r(rr,cc)-pont0));
            hc_half_max = min(hc_half_max1,hc_half_max2);
            
            if hc{rr,cc}<2*hc_half_min
                hc{rr,cc}=2*hc_half_min;
            end

            if hc{rr,cc}>2*hc_half_max
                hc{rr,cc}=2*hc_half_max;
            end
            
            hc_pu{rr,cc} = hc{rr,cc}/(hc_half_max *2);
            
        else
            delta(rr,cc)=(0.5*hc{rr,cc}(1)+sum(hc{rr,cc}(2:end-1))+0.5*hc{rr,cc}(end)-hc_min*(geo.nlay-1))/(Bx0{rr,cc}(1)-Bx0{rr,cc}(end)-hfe_min*(geo.nlay-1)-hc_min*(geo.nlay-1));
            hc_pu{rr,cc}=hc{rr,cc}*(delta(rr,cc)*geo.nlay)/sum(hc{rr,cc});
        end
        % dx (flux carrier design)
        alphad=[0 90-fliplr(geo.alpha)*geo.p 90];                   % 0<=alphad<=90 [° elt]
        r_all=geo.R*xx(rr,cc);
        for ii=1:geo.nlay
            r_all=[r_all Bx0{rr,cc}(ii)+hc{rr,cc}(ii)/2 Bx0{rr,cc}(ii)-hc{rr,cc}(ii)/2];
        end
        r_all=[r_all ArLim(rr,cc)];
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
                B1tmp=ArLim(rr,cc)+hf{rr,cc}(ii);
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
        
        % radial ribs (high speed motors)
        temp.B1k=B1k;
        temp.B2k=B2k;
        temp.Bx0=Bx0;
        temp.yyD2k = [];
        
        geo0.x0 = x0(rr,cc);
        geo0.r = r(rr,cc);
        geo0.RotType = 'Circular';
        geo0.hc = hc{rr,cc};
        if rr==1&&cc==1
            warning('Ribs computed for circular geometry');
        end
        
        [~,geo0] = calc_ribs_rad_Circ(geo0,mat,temp);
        
        pontR(rr,cc) = mean(geo0.pontR);

        % BrPrime for PM-assistance (ribs saturation)
        BrPrime{rr,cc} = Bs*(geo0.pontR/2+geo0.pontT)./(sk{rr,cc}-geo0.pontR/2-geo0.pontT);
        
        BrAvg(rr,cc) = mean(BrPrime{rr,cc});
        BrMax(rr,cc) = max(BrPrime{rr,cc});
        BrMin(rr,cc) = min(BrPrime{rr,cc});

        Abar(rr,cc) = sum(2*sk{rr,cc}.*hc{rr,cc});        
        
        % PM dim (adapt after ribs computation)
        if strcmp(dataSet.TypeOfRotor,'Circular')
            PMdim{rr,cc} = PMdim{rr,cc}-[pontT+geo0.pontR/2;zeros(size(geo0.pontR))];
        elseif strcmp(dataSet.TypeOfRotor,'Seg')|| strcmp(dataSet.TypeOfRotor,'ISeg')
            PMdim{rr,cc} = PMdim{rr,cc}-[geo0.pontR/2;pontT];
        end

        % PM flux linkage
        tmp.nlay    = nlay;
        tmp.Ns      = Ns;
        tmp.g       = g;
        tmp.sk      = sk{rr,cc};
        tmp.hc      = hc{rr,cc};
        tmp.r       = r(rr,cc);
        tmp.dalpha  = geo.dalpha;
        tmp.l       = l;
        tmp.kc      = kc(rr,cc);
        tmp.kw      = kw;
        tmp.PMdim   = PMdim{rr,cc};
        tmp.pontT   = geo0.pontT;
        tmp.pontR   = geo0.pontR;
        tmp.Br      = Br;
        tmp.Bs      = Bs;
        tmp.p       = geo.p;
        tmp.alpha   = alpha;
        tmp.Tfillet = max(geo.RotorFilletTan1,geo.RotorFilletTan2);
        tmp.Rfillet = max(geo.RotorFillet1,geo.RotorFillet2);
                
        fM(rr,cc)   = evalPMfluxSyrmDesign(tmp);
        
 
    end
end

if strcmp(geo.RotType,'SPM')
    % qui devo eliminare tutte le variabili che non servono per SPM (e
    % metterle a NaN) e ricalcolare:
    % - fM
    % - PMdim (per massa PM)
    % - hc (spessore magnete, da formula)
    % - hc_pu
    % - riassegnare dx=dataSet.DepthOfBarrier
    lr    = zeros(m,n);
    PMdim = mat2cell(nan(m,n),m,n);
    hc    = kc*g/1e3./(4/pi*sin(alpha_pu*pi/2)*Br/Bfe./bb-1)*1e3.*kPM;  % PM thickness [mm]
    la    = nan(m,n);
    hcMin = hc;
    skMax = nan(m,n);
    Lfqpu = nan(m,n);
    Lcqpu = nan(m,n);
    pontR = zeros(m,n);
    BrAvg = nan(m,n);
    BrMax = nan(m,n);
    BrMin = nan(m,n);
    Abar  = nan(m,n);
    fM    = kw*Ns*8/pi*sin(alpha_pu*pi/2)*Br./(1+kc*g./(hc)).*R/1e3.*l/1e3/p.*xx;
    for rr=1:m
        for cc=1:n
            hc_pu{rr,cc} = hc(rr,cc)/g;
            dx{rr,cc}    = 1;
        end
    end
end

%%
% d axis

Ht = interp1(mat.Stator.BH(:,1),mat.Stator.BH(:,2),Bfe./kt);
Hy = interp1(mat.Stator.BH(:,1),mat.Stator.BH(:,2),Bfe./ky);
if flag_ks
    %ksat = 1+mu0*pi/2*(Hy*pi/(6*p*q)*(R-ly/2) + Ht*lt + Hy*lr)./(bb*Bfe*kf1.*kc*g);
    % ksat = 1+mu0*(Hy*pi/(6*p*q*n3ph)*(R-ly/2) + Ht*lt + Hy*lr)./(bb*Bfe*kf1.*kc*g);
    ksat = 1+mu0*(Hy*pi/(6*p*q*n3ph)*(R-ly/2) + Ht*lt + Hy*lr)./(bb*Bfe.*kc*g);
else
    ksat=ones(size(xx));
    warning('saturation factor ksat not evaluated!!!')
end

if strcmp(geo.RotType,'SPM')
    id = zeros(m,n);
    Lbase = 6/pi*mu0*(kw*Ns/p).^2.*(R/1e3).*(l/1e3)./(kc.*(g+hc)/1e3).*xx;  % magnetizing unsaturated inductance (same as below, but with increased airgap due to PM thickness)
    Lmd   = Lbase./ksat;    % magnetizing inductance with stator iron saturation partially accounted for
else
    Fmd = 2*(R*1e-3)*(l*1e-3)*kw*Ns*Bfe./p.*xx.*bb;     % flux linkage [Vs]
    id = pi*Bfe*kc*(g*1e-3)*p.*ksat/(mu0*3*kw*Ns).*bb/n3ph;  % id [A]
    %Lmd = Fmd./id/n3ph^2;                                   % d-axis magnetizing inductance [H], with iron saturation
    %Lbase = Lmd.*ksat;  % base inductance for analytical model (unsaturated). The explicit equation is below
    % Lbase = 6/pi*mu0*(kw*Ns/p)^2*(R/1e3)*(l/1e3)./(kc*g/1e3).*xx
    Lbase = 6/pi*mu0*(kw*Ns/p)^2*(R/1e3)*(l/1e3)./(kc*g/1e3).*xx;
    Lmd = Lbase./ksat;
end

% Rated current computation (from kj or J)

per0 = per;
geo0 = geo;

geo0.r = R*xx;
geo0.Aslot = Aslots./(6*p*q*n3ph);
geo0.wt = wt;
geo0.lt = lt;

lend = calc_endTurnLength(geo0);
geo0.lend = lend;

if flag_i0==0
    % compute the rated current from kj (kj=constant)
    per0.J    = NaN;
    per0.kj   = kj;
    per0.Loss = NaN;
    per0.i0   = NaN;
elseif flag_i0==1
    per0.kj   = NaN;
    per0.Loss = NaN;
    per0.J    = J;
    per0.i0   = NaN;
elseif flag_i0==2
    per0.kj   = NaN;
    per0.Loss = NaN;
    per0.J    = NaN;
end

per0 = calc_i0(geo0,per0,mat);
Rs = per0.Rs;
%i0 = real(per0.i0);

kj = per0.kj.*ones(size(xx));
J  = per0.J.*ones(size(xx));
i0 = real(per0.i0).*ones(size(xx));


id(id>loadpu*i0) = loadpu*i0(id>loadpu*i0);

gamma = acos(id./(loadpu*i0));      % current angle from d axis
gamma(gamma<gammaMin) = gammaMin;   % current angle saturation to avoid iq=0
iq    = (loadpu*i0).*sin(gamma);    % q-axis current [A] pk

A = 2*Nbob * loadpu*i0 ./ (r*2*pi/(q*6*p*n3ph)); % linear current density [A/mm] pk

Ad = A.*cosd(gamma);
Aq = A.*sind(gamma);


% copper overtemperature computation
% barriers demagnetization index Ns*i0/hc
dTempCu = zeros(m,n);
NsI0_hc = cell(m,n);
for rr=1:m
    for cc=1:n
        % slot area evaluation
        geo0.r        = geo.R*xx(rr,cc);
        geo0.wt       = wt(rr,cc);
        geo0.lt       = lt(rr,cc);
        geo0.ly       = ly(rr,cc);
        geo0.Aslot    = Aslots(rr,cc)/(6*p*q*n3ph);
        geo0.lend     = lend(rr,cc);
        per0          = per;
        per0.overload = 1;
        per0.Loss     = kj(rr,cc)*(2*pi*R/1000*l/1000);
        
        dTempCu(rr,cc) = temp_est_simpleMod(geo0,per0);
%         % demagnetization index
%         hcTmp = hc{rr,cc};
%         NsI0Tmp = dfQ*Ns*i0(rr,cc);
%         NsI0_hc{rr,cc} = NsI0Tmp./hcTmp;

    end
end
dTempCu = dTempCu-per.temphous;

% tangential ribs effect
if strcmp(geo.RotType,'SPM')
    Lrib = nan(m,n);
    Lmq = Lmd;  % SPM is isotropic
else
    Lrib = 4/pi*kw*Ns*((2*pont+pontR)*1e-3)*(l*1e-3)*Bs./(iq*n3ph);
    Lrib(iq==0) = mean(Lrib(iq~=0));
    Lrib(isnan(Lrib)) = 0;
    Lmq = Lbase.*(Lcqpu+Lfqpu)+Lrib;
end

kdq = 1-Lmq./Lmd;

% T = 3/2*p*(kdq.*Fmd.*iq);

% stator leakage inductance Ls
% [dfs] = staircaseRegular(6*q);                                      % stator staircase
% f = cumsum(dfs);
% sumDf2s = sum(dfs.^2);
d0 = ttd;
d2 = lt- d0 - d1;
beta = c1./c2;
h = (beta.^2-beta.^4/4-log(beta)-0.75)./((1-beta).*(1-beta.^2).^2);
ps = d0./c0 + d1./c0.*1./(c1./c0-1).*log(c1./c0)+d2./c2.*h;

% Lspu = 4/pi*p*kc*g*sumDf2s.*ps./r;                                  % Ls/Lmd: p.u. leakage inductance

% Ls = 2*mu0*(l/1000)*Ns^2/p*ps.*sumDf2s;     % According to Vagati's Tutorial (1994)
Ls = 2*mu0*(l/1000)*Ns^2/p*ps/(q*n3ph);            % According to Lipo's Book (2017)
Lspu = Ls./Lbase;

% power factor
Ld = Lmd+Ls;
Lq = Lmq+Ls;

if strcmp(geo.RotType,'SPM')
    fd = fM;                    % d-axis flux linkage (PM axis)
    fq = n3ph*(Lq.*iq);         % q-axis flux linkage (PM axis)
    ich = fM./Ld;               % characteristic current
    iHWC = abs(fd+j*fq)./Ld;    % Hyper-Worst-Case Peak Short-Circuit current
else
    fd = n3ph*(Ld.*id);         % d-axis flux linkage (SR axis)
    fq = n3ph*(Lq.*iq)-fM;      % q-axis flux linkage (SR axis)
    ich = fM./Lq;               % characteristic current
    iHWC = abs(fd+j*fq)./Lq;    % Hyper-Worst-Case Peak Short-Circuit current
end

T = 3/2*p*n3ph*(fd.*iq-fq.*id);

gamma = atand(iq./id);      % current phase angle from d axis [deg]
delta = atand(fq./fd);      % flux linkage phase angle from d axis [deg]

PF = sind(gamma-delta);     % PF @ gamma (same gamma as torque)

% PFmax = (Ld-Lq)./(Ld+Lq);   % PF @ max PF gamma

% active parts mass computation
mCu = Aslots/1e6.*(l+lend)/1e3*mat.SlotCond.kgm3*kcu;
mPM = nan(size(xx));
if strcmp(mat.LayerMag.MatName,'Air')
    mPM = zeros(size(xx));
else
    if strcmp(geo.RotType,'SPM')
        mPM = pi*((R/1e3*xx).^2-(R/1e3*xx-hc/1e3).^2)*l/1e3*alpha_pu*mat.LayerMag.kgm3;
    else
        for ii=1:numel(xx)
            mPM(ii) = sum(hc{ii}.*sum(PMdim{ii},1))/1e6.*l/1e3*2*(2*p)*mat.LayerMag.kgm3;
        end
    end
end

% Demagnetization indicators
% rev 01: Aqirr from Boazzo
% rev 02: demagnetization volume and minimum PM flux density at i0 and iHWC
tp = pi*r/p;
if Br==0
    Aqirr   = zeros(size(xx));
    Bmin0   = nan(size(xx));
    dPM0    = nan(size(xx));
    BminHWC = nan(size(xx));
    dPMHWC  = nan(size(xx));
else
    if strcmp(geo.RotType,'SPM')
        Bm0pu = nan(m,n);
        Aqirr = nan(m,n);
    else
        % ref: "Design of Ferrite-Assisted Synchronous Reluctance Machines
        % Robust Toward Demagnetization", A. Vagati, 2014
        Bdpu = Bd/Br;
        dalpha = geo.dalpha(end)*pi/180;
        fqn = 1;
        Bm0pu = 1./(1+(skMax/1e3)./(la/1e3).*(g/1e3)./(tp/1e3)*2*pi/dalpha.*sin(dalpha/2));
        Aqirr = pi/4*Br*(la/1e3)./(tp/1e3)/mu0/fqn.*(1-Bdpu./Bm0pu)*2/1e3;                      % multiplied *2 because the full pole is considered here, /1e3 to have it in A/mm
    end

    if dataSet.syrmDesignFlag.demag0
        Bmin0   = Bd*ones(size(xx));
        dPM0    = ones(size(xx));
    else
        Bmin0   = nan(size(xx));
        dPM0    = nan(size(xx));
    end
    if dataSet.syrmDesignFlag.demagHWC
        BminHWC = Bd*ones(size(xx));
        dPMHWC  = ones(size(xx));
    else
        BminHWC = nan(size(xx));
        dPMHWC  = nan(size(xx));
    end
end

if dataSet.syrmDesignFlag.mech==1
    mechStressRad = repmat({repmat(mat.Rotor.sigma_max*10^6,1,nlay)},length(xx(:,1)),length(xx(1,:)));
    mechStressTan = repmat({repmat(mat.Rotor.sigma_max*10^6,1,nlay)},length(xx(:,1)),length(xx(1,:)));
    kmechrad = repmat({repmat(1,1,nlay)},length(xx(:,1)),length(xx(1,:)));
    kmechtan = repmat({repmat(1,1,nlay)},length(xx(:,1)),length(xx(1,:)));
else
    mechStressRad = cell(size(xx));
    mechStressTan = cell(size(xx));
    kmechrad = cell(size(xx));
    kmechtan = cell(size(xx));
end

ktempCuMax    = nan(m,n);
tempCuMax     = dTempCu + per.temphous;
ktempCuMaxAct = nan(m,n);
tempCuMaxAct  = dTempCu + per.temphous;


%% save all the results in the map structure
map.xx            = xx;
map.bb            = bb;
map.T             = T;
map.PF            = PF;
map.id            = id;
map.iq            = iq;
map.fd            = fd;
map.fq            = fq;
map.fM            = fM;
map.wt            = wt;
map.ws            = ws;
map.lt            = lt;
map.la            = la;
map.Ar            = Ar;
map.ArLim         = ArLim;
map.hc_pu         = hc_pu;
map.hc            = hc;
map.dx            = dx;
map.geo           = geo;
map.ly            = ly;
map.ksat          = ksat;
map.Lbase         = Lbase;
map.Lmd           = Lmd;
map.Lcqpu         = Lcqpu;
map.Lfqpu         = Lfqpu;
map.Lspu          = Lspu;
map.Lrpu          = Lrib./Lbase;
map.kdq           = kdq;
map.dTempCu       = dTempCu;
map.Aslots        = Aslots;
map.A             = A;
map.Ad            = Ad;
map.Aq            = Aq;
map.Aqirr         = Aqirr;
map.J             = J;
map.kj            = kj;
map.iAmp          = i0*loadpu;
map.i0            = i0;
map.gamma         = atan2(map.iq,map.id)*180/pi;
map.lend          = lend;
map.BrAvg         = BrAvg;
map.BrMax         = BrMax;
map.BrMin         = BrMin;
map.Abar          = Abar;
map.Rs            = Rs;
map.hcMin         = hcMin;
map.mCu           = mCu;
map.mPM           = mPM;
map.PMdim         = PMdim;
map.NsI0          = map.i0*Ns;
map.F0_Ns         = abs(map.fd+1i*map.fq)/Ns;
map.ich           = ich;
map.iHWC          = iHWC;
map.Bmin0         = Bmin0;
map.dPM0          = dPM0;
map.BminHWC       = BminHWC;
map.dPMHWC        = dPMHWC;
map.mechStressRad = mechStressRad;
map.mechStressTan = mechStressTan;
map.kmechrad      = kmechrad;
map.kmechtan      = kmechtan;
map.NsIch         = Ns.*map.ich;
map.NsIHWC        = Ns.*map.iHWC;
map.kUGO          = fM./(abs(fd+j*fq));
map.ktempCuMax    = ktempCuMax; 
map.tempCuMax     = tempCuMax; 
map.ktempCuMaxAct = ktempCuMaxAct; 
map.tempCuMaxAct  = tempCuMaxAct;
map.NsI0_hc       = NsI0_hc;

map.dataSet       = dataSet;



