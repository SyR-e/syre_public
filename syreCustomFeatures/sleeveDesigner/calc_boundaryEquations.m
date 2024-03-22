function [C_21,p_12,C_22,p_23] =calc_boundaryEquations(par,w,dT_1,dT_2,dT_3,du_12,flag)

%CALC_BOUNDARYEQUATIONS Summary of this function goes here
%   Solve the system of equation in matrix form K*u=G.
%   Every n-th equation can be summarized with the expression 
%   K(n,1)*C_21 + K(n,2)*p_12 + K(n,3)*C_23 + K(n,4)*p_23 = G(n)

% Initialization to zero
K = zeros(4,4);
G = zeros(4,1);

%% Boundary conditions before lift-off

% 1st Equation: K(1,1) * C_21 + p_12 = G(1)
K(1,1) = 1 / par.r_2o;
K(1,2) = 1;
G(1)   = w^2 * (par.rho_2 / 3) * par.r_2o^2;

% 2nd Equation: K(2,1) * C_21 + K(2,4) * p_23 = G(2)
K(2,1) = 1 / par.r_2i;
K(2,4) = 1;
G(2)   = w^2 * (par.rho_2 / 3) * par.r_2i^2;

% 3rd Equation: K(3,1) * C_21 + K(3,2) * p_21 + C_22 = G(3)
K(3,1) = 1/par.E_2 * log(par.r_2o);
K(3,2) = - par.r_1av_loop^2 / (par.E_1 * par.h_1_loop);
K(3,3) = 1;
G(3)   = w^2 * par.rho_1 * (par.r_1av_loop^3 / par.E_1) ....
       + w^2 * (par.rho_2 / (9 * par.E_2)) * par.r_2o^3 ...
       + (par.r_1av_loop * par.a_th1 * dT_1) - (par.r_2o * par.a_th2 * dT_2) ...
       - du_12;

% 4th Equation: K(4,1) * C_21 + C_22 + K(4,4) * p_23 = G(4)
K(4,1) = 1 / par.E_2 * log(par.r_2i);
K(4,3) = 1;
K(4,4) = (1 / par.E_3) * (par.r_3o / (par.r_3o^2 - par.r_3i^2)) * ((1 - par.nu_3) * par.r_3o^2 + (1 + par.nu_3) * par.r_3i^2);
G(4) = w^2 * ((par.rho_2 / (9 * par.E_2)) * par.r_2i^3 + (par.rho_3 / par.E_3) * (3 + par.nu_3) / 8 * par.r_3o * ((par.r_3o^2 + par.r_3i^2) * (1 - par.nu_3) + (1 + par.nu_3) * par.r_3i^2 - ((1 - par.nu_3^2) / (3 + par.nu_3)) * par.r_3o^2)) ...
     - par.r_2i * par.a_th2 * dT_2 + par.r_3o * par.a_th3 * dT_3;

% Solve
sol_BCs = linsolve(K,G);


% assigning values and convert to double

C_21 = double(sol_BCs(1));
p_12 = double(sol_BCs(2));
C_22 = double(sol_BCs(3));
p_23 = double(sol_BCs(4));


%% Boundary conditions after lift-off

if p_23 < 0 && flag == 0

% Initialization to zero
K = zeros(3,3);
G = zeros(3,1);

% 1st Equation: K(1,1) * C_21 + p_12 = G(1)
K(1,1) = 1 / par.r_2o;
K(1,2) = 1;
G(1)   = w^2 * (par.rho_2 / 3) * par.r_2o^2;

% 2nd Equation: K(3,1) * C_21 = G(3)
K(2,1) = 1/par.r_2i;
G(2)   = w^2 * (par.rho_2 / 3) * par.r_2i^2;

% 3rd Equation: K(3,1) * C_21 + K(3,2) * p_21 + C_22 = G(3)
K(3,1) = 1/par.E_2 * log(par.r_2o);
K(3,2) = - par.r_1av_loop^2 / (par.E_1 * par.h_1_loop);
K(3,3) = 1;
G(3)   = w^2 * par.rho_1 * (par.r_1av_loop^3 / par.E_1) ....
       + w^2 * (par.rho_2 / (9 * par.E_2)) * par.r_2o^3 ...
       + (par.r_1av_loop * par.a_th1 * dT_1) - (par.r_2o * par.a_th2 * dT_2) ...
       - du_12;

% Solve
sol_BCs = linsolve(K,G);

% assigning values and convert to double

C_21 = double(sol_BCs(1));
p_12 = double(sol_BCs(2));
C_22 = double(sol_BCs(3));
p_23 = 0;

end
