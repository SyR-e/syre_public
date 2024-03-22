function plot_sleeve_versus_radius(par,in)
%PLOT_SLEEVE_VERSUS_RADIUS Summary of this function goes here
%   Detailed explana.dT_3tion goes here
% plots versus radius for fixed speed, interference and temperature

% conversions
w = 2*pi*(in.n_max/60)/2.5;    % in rad/s. Angular rotor speed.

par.h_1_loop   = par.h_1;
par.r_1o_loop  = par.r_1o;   % in m. Outer sleeve radius.
par.r_1av_loop = par.r_1av;  % in m. Average sleeve radius.

% solve boundary conditions
%    run BCs. 
flag = 0; %In case p_23 < 0, the BCs are changed to loss of contact between inner laminate and PMs + pole piece segment
[C_21,p_12,C_22,p_23] = calc_boundaryEquations(par,w,in.dT_1,in.dT_2,in.dT_3,in.du_12,flag);

% sleeve
s_t1 = p_12 * (par.r_1av / par.h_1) + w^2 * par.rho_1 * par.r_1av^2;    % in Pa. Sleeve circumferential stress.
u_1  = p_12 * (par.r_1av^2 / (par.E_1 * par.h_1)) + w^2 * par.rho_1 * (par.r_1av^3 / par.E_1) + par.r_1av * par.a_th1 * in.dT_1;    % in m. Sleeve radial displacement.

% PM and pole piece segment
r_2 = linspace(par.r_2i, par.r_2o, 1000);   % in m. Radius variable for PM and pole piece segment.

s_r2 = C_21 ./ r_2 - w^2 * (1/3) * par.rho_2 * r_2.^2;    % in Pa. PM and pole piece segment radial stress.
u_2 = (C_21 / par.E_2) .* log(r_2) - w^2 * (par.rho_2 / (9 * par.E_2)) * r_2.^3 + r_2 * par.a_th2 * in.dT_2 + C_22;  % in m. PM and pole piece segment radial displacment.

% inner laminate
r_3 = linspace(par.r_3i, par.r_3o, 1000);   % in m. Radius variable for inner laminate.

s_r3 = (-1) * p_23 * (par.r_3o^2 / (par.r_3o^2 - par.r_3i^2)) * (1 - (par.r_3i ./ r_3).^2) + w^2 * par.rho_3 * ((3 + par.nu_3) / 8) * (par.r_3o^2 + par.r_3i^2 - (par.r_3i * par.r_3o ./ r_3).^2 - r_3.^2);    % in Pa. Inner laminate radial stress.
s_t3 = (-1) * p_23 * (par.r_3o^2 / (par.r_3o^2 - par.r_3i^2)) * (1 + (par.r_3i ./ r_3).^2) + w^2 * par.rho_3 * ((3 + par.nu_3) / 8) * (par.r_3o^2 + par.r_3i^2 + (par.r_3i * par.r_3o ./ r_3).^2 - ((1 + 3*par.nu_3) / (3 + par.nu_3)) * r_3.^2);  % in Pa. Inner laminate circumferential stress.
u_3 = (-1 / par.E_3) * p_23 * (par.r_3o^2 / (par.r_3o^2 - par.r_3i^2)) * ((1 - par.nu_3) .* r_3 + par.r_3i^2 * (1 + par.nu_3) .* (1 ./ r_3)) + w^2 * (par.rho_3 / par.E_3) * ((3 + par.nu_3) / 8) .* r_3 .* ((par.r_3i^2 + par.r_3o^2) * (1 - par.nu_3) + (1 + par.nu_3) .* (par.r_3i .* par.r_3o ./ r_3).^2 - ((1 - par.nu_3^2) / (3 + par.nu_3)) * r_3.^2) + r_3 .* par.a_th3 * in.dT_3;   % in m. Inner laminate radial displacement.


% plotting

% plotting variables and unit conversions
r_1 = linspace(par.r_1i, par.r_1o, 1000);       % in m. Pseudo radius variable for sleeve.
s_t1_plot = s_t1.*linspace(1,1,1000);   % in Pa. Pseudo circumferential stress variable for sleeve.
u_1_plot = u_1.*linspace(1,1,1000);     % in m. Pseudo radial displacement variable for sleeve.

s_t1_plot = s_t1_plot./10^6;  % Pa to MPa
s_r2 = s_r2./10^6;  % Pa to MPa
s_r3 = s_r3./10^6;  % Pa to MPa
s_t3 = s_t3./10^6;  % Pa to MPa
u_1_plot = u_1_plot.*10^6;    % m to um
u_2 = u_2.*10^6;    % m to um
u_3 = u_3.*10^6;    % m to um
r_1 = r_1.*10^3;    % m to mm
r_2 = r_2.*10^3;    % m to mm
r_3 = r_3.*10^3;    % m to mm

% sleeve

figure;
figSetting
subplot(1,2,1)
plot(r_1, s_t1_plot);
xlabel("$r$ in mm");
ylabel("$\sigma_{\rm t1}$ in MPa");
title("Sleeve circumferential stress");
subtitle("$n$ = " + in.n_max/2.5 + " min$^{-1}$, $\Delta u_{12}$ = " + in.du_12*1000 +" mm");
hold off;


subplot(1,2,2)

plot(r_1, u_1_plot);
xlabel("$r$ in mm");
ylabel("$u_{\rm 1}$ in $\mu$m");
title("Sleeve radial displacement");
subtitle("$n$ = " + in.n_max/2.5 + " min$^{-1}$, $\Delta u_{12}$ = " + in.du_12*1000 +" mm");
hold off;

% PM and pole piece segment

figure;
figSetting
subplot(1,2,1)
plot(r_2, s_r2);
xlabel("$r$ in mm");
ylabel("$\sigma_{\rm r2}$ in MPa");
title("PM and pp radial stress");
subtitle("$n$ = " + in.n_max/2.5 + " min$^{-1}$, $\Delta u_{12}$ = " + in.du_12*1000 +" mm");
hold off;


subplot(1,2,2)
plot(r_2, u_2);
xlabel("$r$ in mm");
ylabel("$u_{\rm 2}$ in $\mu$m");
title("PM and pp radial displacement");
subtitle("$n$ = " + in.n_max/2.5 + " min$^{-1}$, $\Delta u_{12}$ = " + in.du_12*1000 +" mm");
hold off;

% inner laminate
figure;
figSetting
subplot(1,3,1)
plot(r_3, s_r3);
title("Inner lam. radial stress");
subtitle("$n$ = " + in.n_max/2.5 + " min$^{-1}$, $\Delta u_{12}$ = " + in.du_12*1000 +" mm");
xlabel("$r$ in mm");
ylabel("$\sigma_{\rm r3}$ in MPa");
hold off;


subplot(1,3,2)
plot(r_3, s_t3);
title("Inner lam. circumferential stress");
subtitle("$n$ = " + in.n_max/2.5 + " min$^{-1}$, $\Delta u_{12}$ = " + in.du_12*1000 +" mm");
xlabel("$r$ in mm");
ylabel("$\sigma_{\rm t3}$ in MPa");
hold off;


subplot(1,3,3)
plot(r_3, u_3);
xlabel("$r$ in mm");
ylabel("$u_{\rm 3}$ in $\mu$m");
title("Inner lam. radial displacement");
subtitle("$n$ = " + in.n_max/2.5 + " min$^{-1}$, $\Delta u_{12}$ = " + in.du_12*1000 +" mm");
hold off;

