function [par,in] = eval_sleeveThickness(dataSet, geo, mat)

fem.res = 0;
fem.res_traf = 0;

% nodes
[rotor,~,geo] = ROTmatr(geo,fem,mat);
[geo,stator,~] = STATmatr(geo,fem);

%general
p   = geo.p;  % Number of pole pairs.

%sleeve
par.h_1   = dataSet.SleeveThickness*10^-3; % in m. sleeve thickness
par.r_1i  = geo.r*10^-3-par.h_1;           % in m. sleeve inner radius
par.r_1o  = par.r_1i + par.h_1;            % in m. Outer sleeve radius.
par.r_1av = (par.r_1i + par.r_1o)/2;       % in m. Average sleeve radius.

%V-shape
a_V  = dataSet.ALPHAdeg*pi/180;
r_Vi = geo.B1k*10^-3; % in m. Radius of inner V-shape gap knee / offset of V-shape gap from rotor center.

%magnet
% PM_points = geo.rotor(geo.rotor(:,9) == 2,:);
h_PM = (geo.B2k-geo.B1k)*cos(a_V)*10^-3;       % in m. Height of buried PM (shorter side without contact with pole piece).
b_PM = dataSet.PMdim(2)*10^-3;   % in m. Width of buried PM (longer side in contact with pole piece).

%pole piece
P0 = [0 0];
P1 = [geo.xxD2k,geo.yyD2k];
P2 = [geo.xxD2k,-geo.yyD2k];
a_pp = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0));
% x_Cpph = 7.7*10^-3;     % in m. Distance from pole piece knee (outer V-shape gap knee) to center of pole piece hole.
r_pph  = 0;             % in m. Radius of pole piece hole.
% r_pph  = 2.4*10^-3;     % in m. Radius of pole piece hole.

%shaft
par.r_3i = dataSet.ShaftRadius*10^-3; % in m. Inner radius of inner laminate (shaft connection).

%% material properties

% sleeve 
par.rho_1  = mat.Sleeve.kgm3;   % in kg/m^3. Sleeve density. 
par.E_1    = mat.Sleeve.E*10^9; % in Pa = N/m^2. Sleeve Young's modulus in fiber direction.
par.nu_1   = 0.28;              % sleeve Poisson number. Typical value from (Zhi, Pingxi, 2019), (Fang, Li, 2019), (Johnson, Hanson, 2021). NOT RELEVANT FOR THE CHOSEN MODELING
par.a_th1  = 0.4*10^-6;         % in 1/K. Sleeve thermal expansion coefficient in fiber direction.
par.s_max1 = 1650*10^6;         % in Pa = N/m^2. Max sleeve allowed stress.
% par.s_max1 = mat.Sleeve.sigma_max*10^6;         % in Pa = N/m^2. Max sleeve allowed stress.

% PMs 
rho_PM = mat.LayerMag.kgm3; % in kg/m^3. PM density. (From SyR-e and BOMATEC datasheet)
E_PM   = 175*10^9;          % in Pa = N/m^2. PM Young's modulus. (Average value between 150 and 200 GPa of datasheet.)
nu_PM  = 0.24;              % PM Poisson number. (DIFFERENT TYPE N38UH, taken from Gao et al 2020.) NOT RELEVANT FOR THE CHOSEN MODELING
a_thPM = 3.5*10^-6;         % in 1/K. PM thermal expansion coefficient in magnetization direction. (Average value between 3 and 4*10^-6/K of datasheet.)

% laminates 
rho_lam = mat.Rotor.kgm3;   % in kg/m^3. Laminates density. (from SyR-e and Gao et al 2020)
E_lam   = mat.Rotor.E*10^9; % in Pa = N/m^2. Laminates Young's modulus. (from SyR-e and Gao et al 2020)
nu_lam  = 0.3;              % Laminates Poisson number. (from Gao et al 2020)
a_thlam = 13*10^-6;         % in 1/K. Laminates thermal expansion coefficient. (from Bicek, Kunc, Zupan 2017)

%% dimensions and geometry

% pole piece 
th_pp = a_pp/2;                                    % in rad. Conversion of angle of V-shape gap.
r_Vo  = geo.B2k*10^-3;                             % in m. Radius of outer V-shape gap knee (pole piece knee).
A_pph = pi*r_pph^2;                                % in m^2. Area of pole piece hole.
h_ppa = par.r_1i*(1 - cos(th_pp));                 % in m. Radial height of arc section of pole piece.
s_ppa = 2*par.r_1i*sin(th_pp);                     % in m. Length of straight face of arc section of pole piece.
h_ppD = par.r_1i - r_Vo - h_ppa;                   % in m. Height of triangle section of pole piece
A_ppD = 0.5*s_ppa*h_ppD;                           % in m^2. Area of triangle section of pole piece.
A_ppa = 0.5*par.r_1i^2*(th_pp - sin(th_pp));       % in m^2. Area of arc section of pole piece.
A_pp  = A_ppD + A_ppa - A_pph;                     % in m^2. Total area of pole piece.
f_h   = A_pph/(A_pp+A_pph);                        % in p.u. Hole fraction of pole piece ("porosity").

% PMs
r_isPMo = sqrt((r_Vi + ...
               (h_PM*tan(a_V) + b_PM)*sin(a_V))^2 + ...
              ((h_PM*tan(a_V) + b_PM)*cos(a_V))^2 );    % in m. Inner surface of PM, outer PM edge radius.
r_isPMi = sqrt((r_Vi + ...
                h_PM*sin(a_V)*tan(a_V))^2 + ...
               (h_PM*sin(a_V))^2);                      % in m. Inner surface of PM, inner PM edge radius.

A_PM    = h_PM*b_PM;                                      % in m^2. Axial cross-section of one PM.
% A_PM = calcAreaShape(PM_points);                       % in m^2. Axial cross-section of one PM.

% combined PMs and pole piece ring segment
par.r_2o = par.r_1i;                      % in m. Outer radius of combined PMs and pole piece ring segment.
par.r_2i = 0.5*(r_isPMo + r_isPMi);       % in m. Inner radius of combined PMs and pole piece ring segment (average PM inner surface radius is taken).
A_2      = th_pp*(par.r_2o^2-par.r_2i^2); % in m^2. Area of combined PMs and pole piece ring segment (with adjusted density).
A_PM_pp  = 2*A_PM + A_pp;                 % in m^2. Area of summarized PMs and pole piece structure.

% inner laminate
par.r_3o = par.r_2i;    % in m. Outer radius of inner laminate ring.

%% material properties

% pole piece
rho_pp = rho_lam;  % in kg/m^3. Pole piece density.
E_ppf  = E_lam;    % in Pa = N/m^2. Pole piece Young's modulus (f = full, without hole influence).
nu_ppf = nu_lam;   % Pole piece Poisson numbe, without hole influence.
a_thpp = a_thlam;  % in 1/K. Pole piece thermal expansion coeffiecent.

% combined PMs and pole piece ring segment
rho_eq    = (rho_PM*2*A_PM + rho_pp*A_pp) / A_PM_pp;                       % in kg/m^3. Equivalent density of summarized PMs and pole piece structure.
par.rho_2 = rho_eq*(A_PM_pp/A_2);                                          % in kg/m^3. Equivalent density of combined PMs and pole piece segment.
f_E       = 1/(1+f_h*3.3775);                                              % in p.u. Lowering factor of Young's modulus with hole influence (randomly distributed "elliptical" shaped holes, Tsukrov, Novak 2001)
E_pp      = f_E*E_ppf;                                                     % in Pa. Young's modulus of pole piece segment with hole influence (randomly distributed "elliptical" shaped holes, Tsukrov, Novak 2001)
par.E_2   = ((2*A_PM/A_PM_pp)*(1/E_PM) + (A_pp/A_PM_pp)*(1/E_pp))^(-1);    % in Pa. Young's modulus of mixed PM and pole piece segment (lower-bound "inverse rule of mixture" for stress perpendicular to "fiber direction").
nu_pp     = (nu_ppf-f_h*(-1.0305))/(1+f_h*3.3775);                         % Poisson number pole piece with hole influence (randomly distributed "elliptical" shaped holes approach, Tsukrov, Novak 2001)
par.nu_2  = ((2*A_PM/A_PM_pp)*(1/nu_PM) + (A_pp/A_PM_pp)*(1/nu_pp))^(-1);  % Poisson number of mixed PM and pole piece segment (lower-bound "inverse rule of mixture" for stress perpendicular to "fiber direction").
par.a_th2 = (2*A_PM/A_PM_pp)*a_thPM + (A_pp/A_PM_pp)*a_thpp;               % in 1/K. Thermal expansion coefficient of mixed PM and pole piece segment without hole influence. (upper-bound "rule of mixture" for more pessimistic thermal expansion)

% inner laminate
par.rho_3 = rho_lam;     % in kg/m^3. Inner laminate density.
par.E_3   = E_lam;       % in Pa = N/m^2. Inner laminate Young's modulus.
par.nu_3  = nu_lam;      % Inner laminate Poisson number.
par.a_th3 = a_thlam;     % in 1/K. Inner laminate thermal expansion coefficient.


%% Plots versus speed for fixed radii, interference and temperature

% inputs
in.n_max      = dataSet.OverSpeed*2.5;              % in min^-1. Maximum rotor speed for plotting.
in.du_12      = dataSet.SleeveInterference*10^(-3); % in m. Assembly interference between sleeve and PMs + pole piece.
in.dT_1       = dataSet.TemperatureRiseAbove20;     % in K. Temperature rise of sleeve above 20 °C.
in.dT_2       = dataSet.TemperatureRiseAbove20;     % in K. Temperature rise of PMs and pole piece above 20 °C.
in.dT_3       = dataSet.TemperatureRiseAbove20;     % in K. Temperature rise of inner laminate above 20 °C.
in.dT         = dataSet.TemperatureRiseAbove20;     % in K. Temperature rise of sleeve above 20 °C.



