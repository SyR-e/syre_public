function [map] = plot_sleeve_design(par,in)
%PLOT_SLEEVE_DESIGN Summary of this function goes here
%   Detailed explanation goes here

%% Speed vectors
% speed variable
n_arr = linspace(0, in.n_max, 70); % in min^-1. Rotor speed.
% unit conversion
w_arr = 2*pi*(n_arr/60);           % in rad/s. Angular rotor speed.

%% Other Inputs
sleeve_max_stress = par.s_max1*10^(-6);          % in MPa. Sleeve max stress.

%% Plot sleeve stress for various sleeve thickness values
h_1_vector = (in.h_1_min:0.1:in.h_1_max)*10^-3;

%Initialize sleeve stress matrix for sleeve thicnkess values
s_t1_h  = zeros(length(h_1_vector), size(w_arr,2));    % in Pa. Sleeve circumferential stress.

for i = 1:length(h_1_vector)
    par.h_1_loop   = h_1_vector(i);
    par.r_1o_loop  = par.r_1i + par.h_1_loop;      % in m. Outer sleeve radius.
    par.r_1av_loop = (par.r_1i + par.r_1o_loop)/2; % in m. Average sleeve radius.

    for k = 1:1:size(w_arr,2)
        w = w_arr(1,k); % indexed angular speed variable

        %run BCs. 
        flag = 0; %In case p_23 < 0, the BCs are changed to loss of contact between inner laminate and PMs + pole piece segment
        [~,p_12,~,~] = calc_boundaryEquations(par,w,0,0,0,0,flag);

        s_t1_h(i,k)  = p_12 * (par.r_1av_loop / par.h_1_loop) + w^2 * par.rho_1 * par.r_1av_loop^2;
    end
end

%Unit convertion
s_t1_h = s_t1_h./10^6;    % Pa to MPa
    

%% Plot sleeve costante

du_12_vector = double((in.du_12_min:0.001:in.du_12_max)*10^-3);

%Initialize
s_t1_p         = zeros(length(du_12_vector), size(w_arr,2));
du_12_break    = zeros(1,length(h_1_vector));

for j = 1:length(h_1_vector)
    par.h_1_loop   = h_1_vector(j);
    par.r_1o_loop  = par.r_1i + par.h_1_loop;      % in m. Outer sleeve radius.
    par.r_1av_loop = (par.r_1i + par.r_1o_loop)/2; % in m. Average sleeve radius.

    for i = 1:length(du_12_vector)
        du_12 = du_12_vector(i);
    
        for k = 1:1:size(w_arr,2)
            w = w_arr(1,k);     % indexed angular speed variable
   
            % run BCs.
            flag = 1;
            [~,p_12,~,~] =calc_boundaryEquations(par,w,in.dT_1,in.dT_1,in.dT_1,du_12,flag);
    
            s_t1_p(i,k,j)  = p_12 * (par.r_1av_loop / par.h_1_loop) + w^2 * par.rho_1 * par.r_1av_loop^2;
        end
    end

    % Calcolo max Sleeve Prestress
    du_12max = 0;
    s_t1_max = 0;
    while s_t1_max <= sleeve_max_stress
        du_12max = du_12max + (sleeve_max_stress - s_t1_max)/sleeve_max_stress*du_12max + 0.0001*10^(-3);
        [~,p_12,~,~] =calc_boundaryEquations(par,0,in.dT_1,in.dT_2,in.dT_3,du_12max,flag);
        s_t1_max  = p_12 * (par.r_1av_loop / par.h_1_loop);
        s_t1_max = s_t1_max/10^6;
    end
    du_12_break(j) = du_12max;

end

%Unit conversion
s_t1_p = s_t1_p./10^6;    % Pa to MPa

%% Calcolo max speed and danger limit

%Initialize to NaN
n_max_sl     = NaN(length(h_1_vector),length(du_12_vector));
n_break      = NaN(length(h_1_vector),length(du_12_vector));
du_12_danger = NaN(length(h_1_vector),length(du_12_vector));

for i = 1:length(h_1_vector)
    for j = 1:length(du_12_vector)
        
        %Lift-off intersection that corrisponds to n_max_sl
        n_max_sl(i,j) = interp1(s_t1_h(i,:) - s_t1_p(j,:,i), n_arr(1,:), 0);
        n_break(i,j)  = interp1(s_t1_p(j,:,i) - sleeve_max_stress, n_arr(1,:), 0);
        
        %Danger area evaluation
        if n_break(i,j) <  n_max_sl(i,j)
            % n_max_sl(i,j) = n_break(i,j);
            du_12_danger(i,j) = du_12_vector(j-1);
        end

        % if s_t1_p(j,:,i) > s_t1_limit(i) 
        %     n_rott(i,j) = n_max_sl(i,j);
        %     n_max_sl(i,j) = NaN;
        % end

        %All the data outside the area of interest become NaN
        if h_1_vector(i) > in.h_1_max*10^(-3)  || du_12_vector(j) > in.du_12_max*10^(-3) 
            n_max_sl(i,j) = NaN;
        end

    end
   du_12_danger_limit(i) = min(du_12_danger(i,:));
end

map.h_1                = h_1_vector;
map.du_12              = du_12_vector;
map.n_max_sl           = n_max_sl - 0.2.*n_max_sl;
map.du_12_break        = du_12_break;
map.du_12_danger_limit = du_12_danger_limit;

    



