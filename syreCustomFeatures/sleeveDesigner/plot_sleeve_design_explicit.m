function plot_sleeve_design_explicit(par,in)
%PLOT_SLEEVE_DESIGN Summary of this function goes here
%   Detailed explanation goes here

% Speed vectors
n_arr = linspace(0, in.n_max, 100); % in min^-1. Rotor speed.
% unit conversion
w_arr = 2*pi*(n_arr/60);         % in rad/s. Angular rotor speed.

% Sleeve max stress
sleeve_max_stress = ones(1,100)*1600;

% Plot settings
figure;
figSetting()
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xlim([0 22500]);
ylim([0 2000])
xlabel("$n$ in min$^{-1}$");
ylabel("$\sigma_{\rm t1}$ in MPa");
title("Sleeve circumferential stress");
jj=0;


%% Plot sleeve thickness costante


h_1_vector = in.h_1_vector/1e3;

%Initialize sleeve stress matrix for sleeve thicnkess values

for ii = 1:length(h_1_vector)
    par.h_1_loop   = h_1_vector(ii);
    par.r_1o_loop  = par.r_1i + par.h_1_loop;      % in m. Outer sleeve radius.
    par.r_1av_loop = (par.r_1i + par.r_1o_loop)/2; % in m. Average sleeve radius.

    for kk = 1:1:size(w_arr,2)
        w = w_arr(1,kk); % indexed angular speed variable

        % run BCs. In case p_23 < 0, the BCs are changed to loss of contact between inner laminate and PMs + pole piece segment
        flag = 0;
        [~,p_12,~,~] = calc_boundaryEquations(par,w,0,0,0,0,flag);
        % [~,p_12,~,~] = calc_boundaryEquations(par,w,in.dT_1,in.dT_2,in.dT_3,in.du_12,flag);

        s_t1_h(ii,kk)  = p_12 * (par.r_1av_loop / par.h_1_loop) + w^2 * par.rho_1 * par.r_1av_loop^2;
    end

    %s_t1_h = s_t1_h(ii,:)./10^6;    % Pa to MPa

    name = sprintf('$h_{sl}$ = %.2g mm',par.h_1_loop*1e3);
    jj = jj+1;
    if par.h_1_loop == par.h_1
        h(jj) = plot(n_arr, s_t1_h(ii,:)/1e6,"LineStyle","-",'DisplayName',name,'LineWidth',2);
    else
        h(jj) = plot(n_arr, s_t1_h(ii,:)/1e6,"LineStyle","-",'DisplayName',name,'LineWidth',1);
    end
    legend(h)
    hold on;

end


%% Plot prestress costante

par.h_1_loop   = par.h_1;
par.r_1o_loop  = par.r_1o;   % in m. Outer sleeve radius.
par.r_1av_loop = par.r_1av;  % in m. Average sleeve radius.

du_12_vector = in.du_12_vector;

for ii = 1:length(du_12_vector)

    du_12 = du_12_vector(ii)*10^-3;

    for kk = 1:1:size(w_arr,2)
        w = w_arr(1,kk);     % indexed angular speed variable

        % run BCs.
        flag = 1;
        [~,p_12,~,~] = calc_boundaryEquations(par,w,in.dT_1,in.dT_2,in.dT_3,du_12,flag);

        s_t1_p(ii,kk)  = p_12 * (par.r_1av_loop / par.h_1_loop) + w^2 * par.rho_1 * par.r_1av_loop^2;
    end

    s_t1_p = s_t1_p(ii,:)./10^6;    % Pa to MPa

    if ii == 1
        jj = jj+1;
        h(jj) = plot(n_arr, s_t1_p,'LineStyle','-',"Color",'k','DisplayName','$du_{12}$','LineWidth',0.5);
    else
        plot(n_arr, s_t1_p,'LineStyle','-',"Color",'k','LineWidth',0.5);
    end

    % Define position to display the text
    ii = round(numel(n_arr)/10);

    % Get the local slope
    d = (s_t1_p(ii+1)-s_t1_p(ii))/(n_arr(ii+1)-n_arr(ii));
    X = diff(get(gca, 'xlim'));
    Y = diff(get(gca, 'ylim'));
    p = pbaspect;
    a = atan(d*p(2)*X/p(1)/Y)*180/pi;

    % Display the text

    label = sprintf('%.2g mm',du_12*10^3);
    text(n_arr(ii), s_t1_p(ii), label, 'BackgroundColor', 'w', 'rotation', a);
    hold on;

end

jj = jj+1;
h(jj) = plot(n_arr, sleeve_max_stress,'LineStyle','--','DisplayName','$\sigma_{max}$','Color','r');

legend (h)


