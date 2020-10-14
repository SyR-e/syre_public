% Copyright 2017
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
%
%    This function is to get Bgap of SPM machine via analytical solution
%    Citation: An improved subdomain model for predicting magnetic field of
%    SPM accounting for tooth-tips,    Z.Q.Zhu,   2011    Trans on Magnetics
%
%    Output: wt: tooth width [mm]
%            ly: yoke length [mm]
%            Bg1: peak fundanmental component of Bgap [T]

function [wt,ly,Bg1] = evalBgapSyrmDesign_SPM(geo,mat,lm_g,x,Bfe)

debug = 0;

p = geo.p;
q = geo.q;
% ps = geo.ps;

if q > 1
    %% DW-SPM
    ps = 1;                                                                 % eval one pole pair
else
    %% CW-SPM
    ps = 2;                                                                 % eval one pole
end

acs = geo.acs;
g = geo.g*1e-3;
R = geo.R*1e-3;
ttd = geo.ttd*1e-3;
Br = mat.LayerMag.Br;
mur = mat.LayerMag.mu;
Beta = geo.betaPMshape;
n3ph = geo.win.n3phase;

mu0 = 4*pi*1e-7;
m = numel(lm_g);
if geo.dx == 1
    PMdir = 'p';
else
    PMdir = 'r';
end

boa = 2*pi/(6*p*q*n3ph)*acs;                                                     % slot opening angle (mech. rad)
bsa = 2*pi/(6*p*q*n3ph)*0.5+0.02;                                                % slot angle (mech. rad), avoid singularity

r = R * x;
lmm = lm_g*g;
lm = mean(lm_g)*g;
Rsb = 0.8*R;
Rr = mean(r) - lm;
Rm = Rr + lm*Beta;
Rs = mean(r) + g;
Rt = Rs + ttd;
Am = geo.phi;
ap = Am/180;

ho = 300;                                                                   % maximum harmonic order
nn = 600*ps;                                                                % number of positions along airgap
nn_PM = round(nn*ap);                                                       % number of PM positions

pc = 2*pi/(6*p*q*n3ph*2);                                                        % half slot span
tau_s = 2*pi/(6*p*q*n3ph);                                                       % slot pitch in radians
tau_PM = ap * pi/p;                                                         % PM span in radians
PM_gap = (1-ap) * pi/p;

rr = Rs - 1/2*g;                                                            % Bgap radius for calculation

kk = 1:1:ho;
mk = 1:1:ho;
nk = 1:1:ho;

%% Initialization
Mrk = zeros(1,ho);
Mak = zeros(1,ho);
ij = find(mod(kk,2*p) == p);                                                 % Eq. 12
%% the magnetization of magnets
if PMdir == 'r'
    %% radial magnetization
    Mrk(1,ij) = 4*p*Br./(ij*pi*mu0).*sin(ij*pi*ap/2/p);                     % radial magnetization k/p = 1,3,5...
else
    %% paralle magnetization
    A1k = sin((kk+1)*ap*pi/(2*p))./((kk+1)*ap*pi/(2*p));                    % Eq. 16
    A2k = sin((kk-1)*ap*pi/(2*p))./((kk-1)*ap*pi/(2*p));                    % Eq. 17
    A2k(1) = 1;
    Mrk(1,ij) = Br/mu0*ap * (A1k(ij) + A2k(ij));                            % Eq. 14  k/p = 1,3,5
    Mak(1,ij) = Br/mu0*ap * (A1k(ij) - A2k(ij));                            % Eq. 15  k/p = 1,3,5
end

En = nk*pi/bsa;                                                             % Eq 32
Fm = mk*pi/boa;                                                             % Eq 36
G1 = diag((Rr/Rm).^kk);                                                     % Eq 23
G2 = diag((Rm/Rs).^kk);                                                     % Eq 54
G3 = diag((Rt/Rsb).^En);                                                    % Eq 33
G4 = diag((Rs/Rt).^Fm);                                                     % Eq 86
Ik = eye(ho);                                                               % Eq 68
K = diag(kk);                                                               % Eq 67
In = eye(ho);
Im = eye(ho);

K11 = Ik + G1.^2;                                                           % Eq 57
K22 = Ik + G1.^2;                                                           % Eq 58
K13 = -G2;                                                                  % Eq 59
K25 = -G2;                                                                  % Eq 60
K14 = -Ik;                                                                  % Eq 61
K26 = -Ik;                                                                  % Eq 62

K31 = Ik - G1.^2;                                                           % Eq 76
K42 = Ik - G1.^2;                                                           % Eq 77
K33 = -mur*G2;                                                              % Eq 78
K45 = -mur*G2;                                                              % Eq 79
K34 = mur*Ik;                                                               % Eq 80
K46 = mur*Ik;                                                               % Eq 81

yy = -2/bsa * En./(Fm'.^2 - En.^2).*(cos(mk'*pi).*sin(En*(bsa + boa)/2) - sin(En*(bsa - boa)/2));

F_m = diag(Fm);                                                             % Eq 99
E_n = diag(En);                                                             % Eq 100
zeta = bsa/boa * yy;                                                        % Eq 115

k97 = yy'*F_m;                                                              % Eq 94
k98 = -yy'*F_m*G4;                                                          % Eq 95
k99 = -E_n*(G3.^2 - In);                                                    % Eq 96

k87 = Im;                                                                   % Eq 112
k88 = G4;                                                                   % Eq 113
k89 = -zeta*(G3.^2 + In);                                                   % Eq 114
K53 = -K;                                                                   % Eq 126
K65 = -K;                                                                   % Eq 126
K54 = G2*K;                                                                 % Eq 127
K66 = G2*K;                                                                 % Eq 127

tm0 = zeros(ho);

%% all slots, considering the rotor rotation
Posts =  linspace(-pi/(2*p)-pc,2*pi-pi/(2*p)-pc,2*6*p*q*n3ph+1);                 % find positions for all slots and teeth center
SlotPos = Posts(2:2:end);                                                   % get all slot center location
ToothPos = Posts(1:2:end);
[~,NearTooth] = min(abs(ToothPos));                                         % find the nearesst tooth to d axis
th_m = ToothPos(NearTooth);                                                 % rotor position, d is on tooth center

%% airgap angular position
theta = linspace(th_m-pi/(2*p),th_m+(2*ps-1)*pi/(p*2),nn);                  % airgap span to eval

%% only slots in rotor theta span are counted
th_s = SlotPos(SlotPos+pc*acs>=theta(1)&SlotPos-pc*acs<=theta(end));
Qs = numel(th_s);                                                           % numbel of slots along theta

%% Magnetization
Mrck = Mrk.*cos(kk*th_m);                                                   % Eq 8
Mrsk = Mrk.*sin(kk*th_m);                                                   % Eq 9
Mack = -Mak.*sin(kk*th_m);                                                  % Eq 10
Mask = Mak.*cos(kk*th_m);                                                   % Eq 11
Y1 = -mu0*((Rr*K*G1 + Rm*Ik)*Mack' - (Rr*G1 + Rm*K)*Mrsk')./(kk.^2-1)';     % Eq 63
Y2 = -mu0*((Rr*K*G1 + Rm*Ik)*Mask' + (Rr*G1 + Rm*K)*Mrck')./(kk.^2-1)';     % Eq 64
Y1(1,1) = mu0*Rm*log(Rm)/2 * (Mack(1)-Mrsk(1));                             % Singularity
Y2(1,1) = mu0*Rm*log(Rm)/2 * (Mask(1)+Mrck(1));                             % Singularity

Y3 = -mu0*(K*(Rm*Ik - Rr*G1)*Mack' - (Rm*Ik - Rr*G1)*Mrsk')./(kk.^2-1)';    % Eq 82
Y4 = -mu0*(K*(Rm*Ik - Rr*G1)*Mask' + (Rm*Ik - Rr*G1)*Mrck')./(kk.^2-1)';    % Eq 83
Y3(1,1) = mu0/2*(Rm*log(Rm) + Rm) * (Mack(1)-Mrsk(1));                      % Singularity
Y4(1,1) = mu0/2*(Rm*log(Rm) + Rm) * (Mask(1)+Mrck(1));                      % Singularity

%% Eq 97
K97 = cell(Qs);
K98 = cell(Qs);
K99 = cell(Qs);
K87 = cell(Qs);
K88 = cell(Qs);
K89 = cell(Qs);
Ftm = cell(Qs);
G4t = cell(Qs);

K97(logical(eye(Qs))) = {k97};
K97(~logical(eye(Qs))) = {zeros(ho)};
K98(logical(eye(Qs))) = {k98};
K98(~logical(eye(Qs))) = {zeros(ho)};
K99(logical(eye(Qs))) = {k99};
K99(~logical(eye(Qs))) = {zeros(ho)};
K87(logical(eye(Qs))) = {k87};
K87(~logical(eye(Qs))) = {zeros(ho)};
K88(logical(eye(Qs))) = {k88};
K88(~logical(eye(Qs))) = {zeros(ho)};
K89(logical(eye(Qs))) = {k89};
K89(~logical(eye(Qs))) = {zeros(ho)};
Ftm(logical(eye(Qs))) = {F_m};
Ftm(~logical(eye(Qs))) = {zeros(ho)};
G4t(logical(eye(Qs))) = {G4};
G4t(~logical(eye(Qs))) = {zeros(ho)};

K97 = cell2mat(K97);
K98 = cell2mat(K98);
K99 = cell2mat(K99);
K87 = cell2mat(K87);
K88 = cell2mat(K88);
K89 = cell2mat(K89);
Ftm = cell2mat(Ftm);
G4t = cell2mat(G4t);

%% Eq 121 & 122
ITA = [];
CSI = [];

for ii = 1:Qs
    ita_si = -1/pi*kk./(Fm'.^2-kk.^2).*(cos(mk'*pi).*sin(kk*th_s(ii) + kk*boa/2)-sin(kk*th_s(ii) - kk*boa/2));
    csi_si =  1/pi*kk./(Fm'.^2-kk.^2).*(cos(mk'*pi).*cos(kk*th_s(ii) + kk*boa/2)-cos(kk*th_s(ii) - kk*boa/2));
    
    [~,temp] = find(Fm'.^2-kk.^2==0);
    ita_si((Fm'.^2-kk.^2)==0) = boa/(2*pi)*cos(kk(temp)*(boa/2-th_s(ii)));
    csi_si((Fm'.^2-kk.^2)==0) = boa/(2*pi)*sin(kk(temp)*(boa/2-th_s(ii)));
    % Eq 132, 133, 134 & 135
    ITA = [ITA; ita_si];
    CSI = [CSI; csi_si];
end
%% C4i, D3i, D4i are independent with slot positions
K57 = ITA'*Ftm*G4t;                                                         % Eq 128
K58 = -ITA'*Ftm;                                                            % Eq 129
K67 = CSI'*Ftm*G4t;                                                         % Eq 130
K68 = -CSI'*Ftm;                                                            % Eq 131

Sigma = 2*pi/boa*ITA;                                                       % Eq 155
Tao = 2*pi/boa*CSI;                                                         % Eq 156

K73 = Sigma;                                                                % Eq 149
K74 = Sigma*G2;                                                             % Eq 150
K75 = Tao;                                                                  % Eq 151
K76 = Tao*G2;                                                               % Eq 152
K77 = -G4t;                                                                 % Eq 153
K78 = -eye(Qs*ho);                                                          % Eq 154

tm1 = zeros(ho,ho*Qs);
tm2 = zeros(ho*Qs,ho);
tm3 = zeros(ho*Qs);

%% Eq. 157
PA = [K11 tm0 K13 K14 tm0 tm0 tm1 tm1 tm1
    tm0 K22 tm0 tm0 K25 K26 tm1 tm1 tm1
    K31 tm0 K33 K34 tm0 tm0 tm1 tm1 tm1
    tm0 K42 tm0 tm0 K45 K46 tm1 tm1 tm1
    tm0 tm0 K53 K54 tm0 tm0 K57 K58 tm1
    tm0 tm0 tm0 tm0 K65 K66 K67 K68 tm1
    tm2 tm2 K73 K74 K75 K76 K77 K78 tm3
    tm2 tm2 tm2 tm2 tm2 tm2 K87 K88 K89
    tm2 tm2 tm2 tm2 tm2 tm2 K97 K98 K99];

tm4 = zeros(size(PA,1)-4*ho,1);
QA = [Y1;Y2;Y3;Y4;tm4];

PAA = sparse(PA);
QAA = sparse(QA);
PB = PAA\QAA;

A2 = PB(2*ho+1:3*ho);
B2 = PB(3*ho+1:4*ho);
C2 = PB(4*ho+1:5*ho);
D2 = PB(5*ho+1:6*ho);

%% calculate flux density, Eq.45 and 46
Brck = -kk.*(A2'/Rs.*(rr/Rs).^(kk-1) + B2'/Rm.*(rr/Rm).^(-kk-1));
Brsk = kk.*(C2'/Rs.*(rr/Rs).^(kk-1) + D2'/Rm.*(rr/Rm).^(-kk-1));
Back = -kk.*(A2'/Rs.*(rr/Rs).^(kk-1) - B2'/Rm.*(rr/Rm).^(-kk-1));
Bask = -kk.*(C2'/Rs.*(rr/Rs).^(kk-1) - D2'/Rm.*(rr/Rm).^(-kk-1));

B2r = full(sum(Brck'.*sin(kk'.*theta) + Brsk'.*cos(kk'.*theta)));           % radial component
B2a = full(sum(Back'.*cos(kk'.*theta) + Bask'.*sin(kk'.*theta)));           % tangential component


%% get Bgap over one pole pair
if ps == 1
    %% DW-SPM
    Bgn = [B2r,-B2r];
    Bgt = [B2a,-B2a];
    theta = [theta,pi/p+theta];
else
    %% CW-SPM
    Bgn = B2r;
    Bgt = B2a;
end


%% for shaped PM
if Beta<1
    xPMco = max(r);
    xPMregular = max(r)-lm + Beta*lm; yPMregular = 0;
    [xPMo,yPMo] = rot_point(xPMregular,yPMregular,ap/2*pi);                 % PM end point
    xArccenter = (xPMco + xPMo - (yPMo.^2./(xPMco-xPMo)))/2;                % find arc center
    Rc = max(r) - xArccenter;
    csi = linspace(-ap*pi/2,ap*pi/2,nn_PM);
    %% calculate PM length along theta
    Lm = (max(r)-Rc)*cos(csi) + sqrt(Rc^2-(max(r)*sin(csi)-Rc*sin(csi)).^2)-max(r)+lm;
    gg = lm + g - Lm;                                                       % airgap length
    
    AirPortionShape_1 = Bgn(1:round((nn - nn_PM)/2/ps));
    AirPortionShape_2 = Bgn(round((nn + nn_PM)/2/ps)+1:round((3*nn - nn_PM)/2/ps));
    AirPortionShape_3 = Bgn(round((3*nn + nn_PM)/2/ps)+1:end);
    PMPortionShape_1 = Bgn(round((nn - nn_PM)/2/ps)+1:round((nn + nn_PM)/ps/2));
    PMPortionShape_2 = Bgn(round((3*nn - nn_PM)/2/ps)+1:round((3*nn + nn_PM)/ps/2));
    
    ratio_shaped = Lm*(Lm(1)+gg(1)*mur)./((Lm+gg*mur)*Lm(1));
    PMPortionShape1 = PMPortionShape_1.*ratio_shaped;
    PMPortionShape2 = PMPortionShape_2.*ratio_shaped;
    Bgn = [AirPortionShape_1,PMPortionShape1,AirPortionShape_2,PMPortionShape2,AirPortionShape_3];
end

%% find PM portion and apply to other lm/g over one pole pair
AirPortion_1 = Bgn(1:round((nn - nn_PM)/ps/2)).*ones(m,1);
AirPortion_2 = Bgn(round((nn + nn_PM)/ps/2)+1:round((3*nn - nn_PM)/ps/2)).*ones(m,1);
AirPortion_3 = Bgn(round((3*nn + nn_PM)/ps/2)+1:end).*ones(m,1);
PMPortion_1 = Bgn(round((nn - nn_PM)/ps/2)+1:round((nn + nn_PM)/ps/2));
% PMPortion_1(PMPortion_1>PMPortion_1(end/2)) = PMPortion_1(end/2);
PMPortion_2 = Bgn(round((3*nn - nn_PM)/ps/2)+1:round((3*nn + nn_PM)/ps/2));
ratio = lmm.*(lm+g*mur)./((lmm+g*mur)*lm);
PMPortion1 = PMPortion_1.*ratio';
PMPortion2 = PMPortion_2.*ratio';
%% one pole pair distribution
Bg_pp = [AirPortion_1,PMPortion1,AirPortion_2,PMPortion2,AirPortion_3];
%% ps = 1
% AirPortion_1 = Bgn(1:round((nn - nn_PM)/ps/2)).*ones(m,1);
% AirPortion_2 = Bgn(round((nn + nn_PM)/ps/2)+1:end).*ones(m,1);
% PMPortion_1 = Bgn(round((nn - nn_PM)/ps/2)+1:round((nn + nn_PM)/ps/2));
% ratio = lmm.*(lm+g*mur)./((lmm+g*mur)*lm);
% PMPortion1 = PMPortion_1.*ratio';
% Bg_pole = [AirPortion_1,PMPortion1,AirPortion_2];
%% get average Bgap, only PM portion is counted
Bg_avg = mean(abs(PMPortion1),2);

%% five different cases
if q > 1 || q == 4/5 || q==4/7
        %% DW-SPM
    % find most loaded tooth position
    [~,mt] = min(abs(th_s-th_m));
    Bg_tooth= Bg_pp(:,theta>=th_s(mt)&theta<=th_s(mt+1));
    B_pm = mean(Bg_tooth,2);
    
    %% sizing
    wt = 2*pi*r/(6*p*q*n3ph).*B_pm / Bfe;
    ly = pi*r.*Bg_avg*ap/(2*p)/Bfe;  
else    
    %% CW-SPM
    if tau_PM < tau_s
        
        if PM_gap <  tau_s * acs/2
            %% when two opposite PMs are so close, only PM part against tooth is counted
            tm7 = tau_s*(1-acs)/2+th_m;
            tm8 = -tau_s*(1-acs)/2+th_m;
            Bg_tooth= Bg_pp(:,theta>=tm8&theta<=tm7);
            B_pm = mean(Bg_tooth,2);
            
            %% sizing
            wt = B_pm/Bfe*tau_s*(1-acs)*r;
        else
            %% PM area < slot pitch
            B_pm = Bg_avg;
            
            %% sizing
            wt = B_pm/Bfe*tau_PM*r;
        end
%     elseif  tau_PM >= tau_s && tau_PM < tau_s*(1+acs)
    else
        %% PM area > one slot span & less than two slots span
        Bgt_pole = Bgt(1:end/2);
        Bg_pole = Bg_pp(:,1:end/2);
        
        [~,maxBgt_p] = max(Bgt_pole(1:end/2));                              % find max Bgt
        
        %% find zero point next to max Bgt point
        tm5 = sign(Bgt(1:end/2).*[Bgt(2:end/2),Bgt(1)]);
        [~,tm6] = find (tm5 < 0);
        tm7 = tm6(sum(tm6-maxBgt_p<0));
        tm8 = size(Bgt_pole,2) - tm7;
        Bg_tooth= Bg_pole(:,tm7:tm8+1);
        Active_PM_ratio = (tm8+1-tm7)/size(Bg_pole,2);
        B_pm = mean(Bg_tooth,2);
        Active_PM_pitch = Active_PM_ratio * pi*(r+g)/p;
        %% sizing
        wt = B_pm/Bfe*Active_PM_pitch;
    end
    %% sizing
    ly = wt/2;
end

wt = wt*1e3;
ly = ly*1e3;

%% calculate Bg1
LL = size(Bg_pp,2);
Bg1 = size(1,m);
for mm = 1:m
    Y = fft(Bg_pp(mm,:));
    P2 = abs(Y/LL);
    Bg1(mm) = 2*P2(2);                                                      % get Bg1 from fft of Bg_pole
end

if (debug)
    figure(11)
    hold on
    grid on
    plot(theta,Bgn,'r');
    plot(theta,Bgt,'b')
    xlim([theta(1),theta(end)])
    keyboard
end
