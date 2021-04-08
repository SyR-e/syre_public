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

function [map] = syrmDesign_SPM(dataSet)
% 
% [map] = syrmDesign_SPM(dataSet)
% 
% Preliminary design for SPM machines. References:
% - ECCE 2017 paper ("Parametric design method for SPM machines including rounded PM shape")

%% inputs
mu0 = 4e-7*pi;                                                              % air permeability

Bfe = dataSet.Bfe;                                                          % steel loading (yoke flux density [T])
kt = dataSet.kt;


[~, ~, geo, per, mat] = data0(dataSet);

R = geo.R;
p = geo.p;
q = geo.q;
acs = geo.acs;
g = geo.g;
l = geo.l;
avv = geo.win.avv;
kcu = geo.win.kcu;
Ns = geo.win.Ns;
Nbob = geo.win.Nbob;
kracc = geo.win.kracc;
Qs = geo.Qs;
phi = geo.phi;
% ns = geo.ns;
ttd = geo.ttd;
tta = geo.tta;
RaccordoFC = geo.SFR;
nlay = geo.nlay;
pont0 = geo.pont0;
pont = geo.pont0;%+min(geo.pont);
n3ph = geo.win.n3phase;



Br = mat.LayerMag.Br;
mur = mat.LayerMag.mu;
Loss = per.Loss;
tempcu = per.tempcu;

if isfield(mat.LayerMag,'temp')
    tempPM = per.tempPP;
    BrTemp = mat.LayerMag.temp;
    Br = interp1(BrTemp.temp,BrTemp.Br,tempPM);
else
    Br = mat.LayerMag.Br;
end


%flags for SyRM design
flag_sd = 0;    % flag_sd=1 --> use subdomain model to compute Fmd,wt,ly (Chao Lu PhD dissertation)
                % flag_sd=0 --> use pure analytical equation to compute Fmd,wt,ly (ECCE 2017)

% design domain according to x and lm
m = 31; n = 21;                                                     % m x n grid of evaluated machines
lm_g = linspace(dataSet.bRange(1),dataSet.bRange(2),m);
x = linspace(dataSet.xRange(1),dataSet.xRange(2),n);                % rotor/stator split factor
ap = geo.phi/180;                                                   % PM/pole ratio

% parametric analysis: design domain (x,b)
[xx,lm_gp] = meshgrid(x,lm_g);

beta = geo.betaPMshape;
r = R * xx;
rocu = 17.8*(234.5 + tempcu)/(234.5+20)*1e-9;                    % resistivity of copper [Ohm m]
ssp = (r+g) * pi/(3*p*q);                                           % stator slot pitch (x)
sso = ssp * acs;                                                    % stator slot opening (x)
Ar = r*(sqrt(2)-1);                                                 % shaft radius (x) [mm]
[kw, ~] = calcKwTh0(geo);
kw = kw(1);

if ~flag_sd
    %Bg calculation
    lm = lm_gp*g;
    xPMco = r;
    xPMregular = r-lm + beta*lm;
    yPMregular = 0;
    [xPMo,yPMo] = rot_point(xPMregular,yPMregular,phi/2*pi/180);        % PM edge point
    xArccenter = (xPMco + xPMo - (yPMo.^2./(xPMco-xPMo)))/2;            % find arc center location
    Rc = r - xArccenter;

    csi = linspace(-phi*pi/180/2,phi*pi/180/2,300);
    air = zeros(size(csi,1),ceil((180*size(csi,2)/phi-size(csi,2))/2)); % the size of no mag zone relates to Am
    for mm = 1:m
        for nn = 1:n
            Lm{mm,nn} = (r(mm,nn)-Rc(mm,nn))*cos(csi) + sqrt(Rc(mm,nn)^2-(r(mm,nn)*sin(csi)-Rc(mm,nn)*sin(csi)).^2)-r(mm,nn)+lm(mm,nn);
            G{mm,nn} = lm(mm,nn) +geo.g - Lm{mm,nn};
            kc{mm,nn} = ssp(mm,nn)./(ssp(mm,nn)-2/pi*G{mm,nn}.*(sso(mm,nn)./G{mm,nn}.*atan(sso(mm,nn)./(2*G{mm,nn}))-log(1+(sso(mm,nn)./(2*G{mm,nn})).^2)));
            Bg{mm,nn} = Lm{mm,nn}./G{mm,nn}./(Lm{mm,nn}./G{mm,nn}+kc{mm,nn}*mur)*Br;
            temp{mm,nn} = [air,Bg{mm,nn},air];
%             Bg_avg(mm,nn) = mean(Bg{mm,nn});
            Bg_avg(mm,nn) = mean(temp{mm,nn});
            Bg_pole{mm,nn} = [temp{mm,nn},-temp{mm,nn}];
            L = length(Bg_pole{mm,nn});
            Y = fft(Bg_pole{mm,nn});
            P2 = abs(Y/L);
            Bg1(mm,nn) = 2*P2(2);                                       % get Bg1 from fft of Bg
            Bg_max(mm,nn) = max(Bg{mm,nn});
        end
    end
    %wt = 2*pi*r.*Bt_max/(6*p*q)/Bfe;
    %ly = pi*r.*Bg_avg*ap/(2*p)/Bfe;                                     % Bianchi 'Theory and design of fractional-slot pm machines'(7.1)
    if q < 1            
        ly = pi*r.*Bg_avg*phi/180/(2*p)/Bfe*(1-acs);                    
        wt = 2*pi*r.*Bg_avg/(6*p*q*n3ph)/Bfe;                                % Hanselman 'Brushless PM machine design' (9.4)
    else
%         ly = pi*r.*Bg_avg*phi/180/(2*p)/Bfe;                            % Hanselman 'Brushless PM machine design' (9.7)
        ly = R/p*pi/2*Bg_avg/Bfe.*xx;
%         wt = 2*pi*r.*Bg1/(6*p*q)/Bfe;                                % Hanselman 'Brushless PM machine design' (9.4)
    wt = 2*pi*R/(6*p*q*n3ph)*Bg_max/Bfe.*xx.*kt;    
    end
else    %sub-domain model (Chao Lu)
    [wt,ly,Bg1] = evalBgapSyrmDesign_SPM(geo,mat,lm_g,x,Bfe);
end

if Br ==0
    h = errordlg('Please use a real magnet material and define Br in Other Options tab');
    uiwait(h);
    return
end


%  
lt = R - r -g - ly;                                                 % slot length (x,lm) [mm]
lt(lt<2*ttd)=NaN;

% d axis
if flag_sd
    Fmd = 2*r.*Bg1'*l*Ns*kw/p*1e-6;
else
    Fmd = pi*l*Ns/(sqrt(3)*p)* Bg1.*r*1e-6;
end
% stator design
if q<1
    lend = 0.5*(wt+pi*(r+lt/2)*sin(pi/(6*p*q*n3ph)));
else
    lend = 2*lt+(r+g+lt/2)*2*pi/(p*q*n3ph);                              % end-turn length [mm]
end

%% calculate slot area (regualr region subtract redundant region of fillet radius)
Aslots  = zeros(m,n);
dTempCu = zeros(m,n);
for ii=1:m
    for jj=1:n
        % slot area evaluation
        geo0 = geo;
        geo0.r  = r(ii,jj);
        geo0.lt = lt(ii,jj);
        geo0.wt = wt(ii,jj);
        try
            [tmp1,~] = drawSlot(geo0);
            Aslots(ii,jj) = 6*geo.p*geo.q*tmp1.Aslot;
        catch
            Aslots(ii,jj) = NaN;
        end
        
        % copper overtemperature
        geo0=geo;
        geo0.r  = geo.R*xx(ii,jj);
        geo0.wt = wt(ii,jj);
        geo0.lt = lt(ii,jj);
        geo0.ly = ly(ii,jj);
        
        dTempCu(ii,jj) = temp_est_simpleMod(geo0,per);
    end
end

dTempCu=dTempCu-per.temphous;

% Aslots = 2 * area_half_slot *6*p*q;
Aslots(Aslots<0)=NaN;

kj = Loss/(2*pi*R*l)*1e6;                                           % specific loss (x,lm) [W/m2]
K = sqrt(kcu*kj/rocu*l./(l+lend));                                  % factor K [] (x,lm)
i0 = pi/(3*Ns)*(R/1000)^1.5*K.*sqrt(Aslots/(pi*R^2));               % rated current i0 [A] pk
loadpu = dataSet.CurrLoPP;                                          % current load in p.u. of i0

Nbob  = Ns/p/q/2;                                                   % conductors in slot per layer
J = 2*Nbob * loadpu*i0 ./ (Aslots/(q*6*p)*kcu);                     % current density in copper [A/mm2] pk
A = 2*Nbob * loadpu*i0 ./ (r*2*pi/(q*6*p));                         % linear current density [A/mm] pk

%% Inductance calculation
kc = ssp./(ssp-2/pi*g*(sso/g.*atan(sso/(2*g))-log(1+(sso/(2*g)).^2)));   % Carter coefficient (x)

% megnetization inductance
% use the formula 3.110 of Pyrhonen: total Lm
Lmd = (6/pi*mu0)/p^2*r.*l./(lm_gp*g+g*kc).*(Ns*kw)^2;               %[mH] Juha Pyrhonen 'Design of rotating electrical machines' (3.110)

% slot leakage inductance, dependent on slot shape
h1 = ttd;
b1 = (r+g)*pi/(3*p*q)*acs;
h2 = 0;
h4 = lt-h1-h2;
b4 = pi*(r+g+0.5*lt)./(3*p*q)-wt;
Lslot_self = 12*(h4./b4/3+h1./b1)/(6*p*q)*mu0*l*(Ns*kw)^2;          %[mH] Juha Pyrhonen ' Design of rotating electrical machines' (4.30)
%% mutual inductance is included
Lslot_mutual = 2*(h4./b4/3+h1./b1)/(6*p*q)*mu0*l*(Ns*kw)^2;
Lslot = Lslot_self + Lslot_mutual;                                  %[mH] El-Refaie "Winding Inductances of Fractional Slot Surface-Mounted Permanent Magnet Brushless Machines," (17)

% tip leakage inductance, short pitching is neglected
Ltip = 5*((g+lm_gp*g/mur)./b1)./(5+4*(g+lm_gp*g/mur)./b1)*3*4/(6*p*q)*mu0*l*(Ns*kw)^2;  % [mH] Juha Pyrhonen ' Design of rotating electrical machines' (4.63)
iq = i0*loadpu;                                                     % q-axis current [A] pk
T = 3/2*p*Fmd.*iq;

%% PF calculation
Ld = (Lmd + Lslot + Ltip)*1e-3;
Lq = Ld;
PF = Fmd./sqrt(Fmd.^2 + (Lq.* iq).^2);


%% save all the results in the map structure
map.xx      = xx;
map.bb      = lm_gp;
map.T       = T;
map.PF      = PF;
map.id      = zeros(m,n);
map.iq      = iq;
map.fd      = Fmd;
map.fM      = Fmd;
map.fq      = Lq.*iq;
map.wt      = wt;
map.lt      = lt;
map.Ar      = Ar;
map.geo     = geo;
map.ly      = ly;
map.Lmd     = Lmd;
map.Ltip    = Ltip;
map.Lslot   = Lslot;
map.Aslots  = Aslots;
map.iAmp    = i0*loadpu;
map.gamma   = atan2(map.iq,map.id)*180/pi;
map.A       = A;
map.J       = J;
map.dtempCu = dTempCu;
map.dataSet = dataSet;



