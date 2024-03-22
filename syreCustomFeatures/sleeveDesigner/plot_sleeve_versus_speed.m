function plot_sleeve_versus_speed(par,in)
%PLOT_SLEEVE_VERSUS_SPEED Summary of this function goes here
%   Detailed explanation goes here

% plots versus speed for fixed radii, interference and temperature

% speed variable
n_arr = linspace(0, in.n_max, 200);    % in min^-1. Rotor speed.
% unit conversion
w_arr = 2*pi*(n_arr/60);    % in rad/s. Angular rotor speed.

% generate solution arrays
s_t1  = zeros(1, size(w_arr,2));    % in Pa. Sleeve circumferential stress.
u_1   = zeros(1, size(w_arr,2));    % in m. Sleeve radial displacement.
s_r2o = zeros(1, size(w_arr,2));    % in Pa. PM + pole piece segment r_1avradial stress at its outer radius.
s_r2i = zeros(1, size(w_arr,2));    % in Pa. PM + pole piece segment radial stress at its inner radius.
u_2o  = zeros(1, size(w_arr,2));    % in m. PM + pole piece segment radial displacement at its outer radius.
u_2i  = zeros(1, size(w_arr,2));    % in m. PM + pole piece segment radial displacement at its inner radius.
s_r3o = zeros(1, size(w_arr,2));    % in Pa. Inner laminate radial stress at its outer radius.
s_r3i = zeros(1, size(w_arr,2));    % in Pa. Inner laminate radial stress at its inner radius.
s_t3o = zeros(1, size(w_arr,2));    % in Pa. Inner laminate circumferential stress at its outer radius.
s_t3i = zeros(1, size(w_arr,2));    % in Pa. Inner laminate circumferential stress at its inner radius.
u_3o  = zeros(1, size(w_arr,2));    % in m. Inner laminate radial displacement at its outer radius.
u_3i  = zeros(1, size(w_arr,2));    % in m. Inner laminate radial displacement at its inner radius.

for k = 1:1:size(w_arr,2)
    w = w_arr(1,k);     % indexed angular speed variable

    par.h_1_loop   = par.h_1;
    par.r_1o_loop  = par.r_1o;   % in m. Outer sleeve radius.
    par.r_1av_loop = par.r_1av;  % in m. Average sleeve radius.


    % run BCs. 
    flag = 0; %In case p_23 < 0, the BCs are changed to loss of contact between inner laminate and PMs + pole piece segment
    [C_21,p_12,C_22,p_23] = calc_boundaryEquations(par,w,in.dT_1,in.dT_2,in.dT_3,in.du_12,flag);
    
    s_t1(1,k)  = p_12 * (par.r_1av / par.h_1) + w^2 * par.rho_1 * par.r_1av^2;
    u_1(1,k)   = p_12 * (par.r_1av^2 / (par.E_1 * par.h_1)) + w^2 * par.rho_1 * (par.r_1av^3 / par.E_1) + par.r_1av * par.a_th1 * in.dT_1;
    s_r2o(1,k) = C_21 / par.r_2o - w^2 * (1/3) * par.rho_2 * par.r_2o^2;
    s_r2i(1,k) = C_21 / par.r_2i - w^2 * (1/3) * par.rho_2 * par.r_2i^2;
    u_2o(1,k)  = (C_21 / par.E_2) * log(par.r_2o) - w^2 * (par.rho_2 / (9 * par.E_2)) * par.r_2o^3 + par.r_2o * par.a_th2 * in.dT_2 + C_22;
    u_2i(1,k)  = (C_21 / par.E_2) * log(par.r_2i) - w^2 * (par.rho_2 / (9 * par.E_2)) * par.r_2i^3 + par.r_2i * par.a_th2 * in.dT_2 + C_22;
    s_r3o(1,k) = (-1) * p_23 * (par.r_3o^2 / (par.r_3o^2 - par.r_3i^2)) * (1 - (par.r_3i / par.r_3o)^2) + w^2 * par.rho_3 * ((3 + par.nu_3) / 8) * (par.r_3o^2 + par.r_3i^2 - (par.r_3i * par.r_3o / par.r_3o)^2 - par.r_3o^2);
    s_r3i(1,k) = (-1) * p_23 * (par.r_3o^2 / (par.r_3o^2 - par.r_3i^2)) * (1 - (par.r_3i / par.r_3i)^2) + w^2 * par.rho_3 * ((3 + par.nu_3) / 8) * (par.r_3o^2 + par.r_3i^2 - (par.r_3i * par.r_3o / par.r_3i)^2 - par.r_3i^2);
    s_t3o(1,k) = (-1) * p_23 * (par.r_3o^2 / (par.r_3o^2 - par.r_3i^2)) * (1 + (par.r_3i / par.r_3o)^2) + w^2 * par.rho_3 * ((3 + par.nu_3) / 8) * (par.r_3o^2 + par.r_3i^2 + (par.r_3i * par.r_3o / par.r_3o)^2 - ((1 + 3*par.nu_3) / (3 + par.nu_3)) * par.r_3o^2);
    s_t3i(1,k) = (-1) * p_23 * (par.r_3o^2 / (par.r_3o^2 - par.r_3i^2)) * (1 + (par.r_3i / par.r_3i)^2) + w^2 * par.rho_3 * ((3 + par.nu_3) / 8) * (par.r_3o^2 + par.r_3i^2 + (par.r_3i * par.r_3o / par.r_3i)^2 - ((1 + 3*par.nu_3) / (3 + par.nu_3)) * par.r_3i^2);
    u_3o(1,k)  = (-1 / par.E_3) * p_23 * (par.r_3o^2 / (par.r_3o^2 - par.r_3i^2)) * ((1 - par.nu_3) * par.r_3o + par.r_3i^2 * (1 + par.nu_3) * (1 / par.r_3o)) + w^2 * (par.rho_3 / par.E_3) * ((3 + par.nu_3) / 8) * par.r_3o * ((par.r_3i^2 + par.r_3o^2) * (1 - par.nu_3) + (1 + par.nu_3) * (par.r_3i * par.r_3o / par.r_3o)^2 - ((1 - par.nu_3^2) / (3 + par.nu_3)) * par.r_3o^2) + par.r_3o * par.a_th3 * in.dT_3;
    u_3i(1,k)  = (-1 / par.E_3) * p_23 * (par.r_3o^2 / (par.r_3o^2 - par.r_3i^2)) * ((1 - par.nu_3) * par.r_3i + par.r_3i^2 * (1 + par.nu_3) * (1 / par.r_3i)) + w^2 * (par.rho_3 / par.E_3) * ((3 + par.nu_3) / 8) * par.r_3i * ((par.r_3i^2 + par.r_3o^2) * (1 - par.nu_3) + (1 + par.nu_3) * (par.r_3i * par.r_3o / par.r_3i)^2 - ((1 - par.nu_3^2) / (3 + par.nu_3)) * par.r_3i^2) + par.r_3i * par.a_th3 * in.dT_3;

end

% plotting

% unit conversions
s_t1 = s_t1./10^6;    % Pa to MPa
s_r2o = s_r2o./10^6;  % Pa to MPa
s_r2i = s_r2i./10^6;  % Pa to MPa
s_r3o = s_r3o./10^6;  % Pa to MPa
s_r3i = s_r3i./10^6;  % Pa to MPa
s_t3o = s_t3o./10^6;  % Pa to MPa
s_t3i = s_t3i./10^6;  % Pa to MPa
u_1 = u_1.*10^6;      % m to um
u_2o = u_2o.*10^6;    % m to um
u_2i = u_2i.*10^6;    % m to um
u_3o = u_3o.*10^6;    % m to um
u_3i = u_3i.*10^6;    % m to um

% sleeve
figure;
subplot(1,2,1)
figSetting
plot(n_arr, s_t1);
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xlim([0 22500]);
xlabel("$n$ in min$^{-1}$");
ylabel("$\sigma_{\rm t1}$ in MPa");
title("Sleeve circumferential stress, $\Delta u_{12}$ = " + in.du_12*1000 +" mm");
hold off;


subplot(1,2,2)
plot(n_arr, u_1);
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xlim([0 22500]);
xlabel("$n$ in min$^{-1}$");
ylabel("$u_{\rm 1}$ in $\mu$m");
title("Sleeve radial displacement, $\Delta u_{12}$ = " + in.du_12*1000 +" mm");
hold off;

% PM and pole piece segment
figure;
figSetting
subplot(1,2,1)
plot(n_arr, s_r2o);
plot(n_arr, s_r2i);
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xlim([0 22500]);
xlabel("$n$ in min$^{-1}$");
ylabel("$\sigma_{\rm r2}$ in MPa");
legend("$\sigma_{\rm r2}(r_{\rm 2o})$", "$\sigma_{\rm r2}(r_{\rm 2i})$");
title("PM and pp radial stress, $\Delta u_{12}$ = " + in.du_12*1000 +" mm");
hold off;

subplot(1,2,2)
plot(n_arr, u_2o);
plot(n_arr, u_2i);
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xlim([0 22500]);
xlabel("$n$ in min$^{-1}$");
ylabel("$u_{\rm 2}$ in $\mu$m");
legend("$u_{\rm 2}(r_{\rm 2o})$", "$u_{\rm 2}(r_{\rm 2i})$");
title("PM and pp radial displacement, $\Delta u_{12}$ = " + in.du_12*1000 +" mm");
hold off;

% inner laminate
figure;
figSetting
subplot(1,3,1)
plot(n_arr, s_r3o);
plot(n_arr, s_r3i);
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xlim([0 22500]);
title("Inner lam. radial stress, $\Delta u_{12}$ = " + in.du_12*1000 +" mm");
xlabel("$n$ in min$^{-1}$");
ylabel("$\sigma_{\rm r3}$ in MPa");
legend("$\sigma_{\rm r3}(r_{\rm 3o})$", "$\sigma_{\rm r3}(r_{\rm 3i})$");
hold off;

subplot(1,3,2)
plot(n_arr, s_t3o);
plot(n_arr, s_t3i);
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xlim([0 22500]);
title("Inner lam. circumferential stress, $\Delta u_{12}$ = " + in.du_12*1000 +" mm");
xlabel("$n$ in min$^{-1}$");
ylabel("$\sigma_{\rm t3}$ in MPa");
legend("$\sigma_{\rm t3}(r_{\rm 3o})$", "$\sigma_{\rm t3}(r_{\rm 3i})$");
hold off;

subplot(1,3,3)
plot(n_arr, u_3o);
plot(n_arr, u_3i);
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xlim([0 22500]);
xlabel("$n$ in min$^{-1}$");
ylabel("$u_{\rm 3}$ in $\mu$m");
legend("$u_{\rm 3}(r_{\rm 3o})$", "$u_{\rm 3}(r_{\rm 3i})$");
title("Inner lam. radial displacement, $\Delta u_{12}$ = " + in.du_12*1000 +" mm");
hold off;

