%% VERSIONE 20 11 2011
% FluxElab
% PostProcess Flux_UVW --> Flux_dq

% phase flux
% time sequence is POSITIVE
% physical sequence is POSITIVE

% coil positive is towards airgap 


% 
% % eliminates time = 0
% Fluxa = Flux_1sim(:,1);
% Fluxb = Flux_1sim(:,2);
% Fluxc = Flux_1sim(:,3);

% rotor mechanical position (RunS_SetParCase_new)
% t = 0, fase corrente in base a gamma, posizione rotore 0
pos_0 = 0;
pos_Mov = time/1000 * per.EvalSpeed * pi/30 + pos_0;  % rotor position (rad mec)
pos_Mov = pos_Mov(1:end)';

%% posizione del rotore (rad elt)
if exist('th0')
    % casi successivi feb 2010 (+90 perchè Waveform in Magnet\coil è un sin)
    theta = th0(1) * pi/180 + pos_Mov * geo.p;
else
    % casi precedenti feb 2010
    % d_axis Vs alpha_axis mechanical position
    theta = pi/2 + pos_Mov * geo.p;  % mec rad
end
% dq flux
% FluxAlpha = 2/3 * (Fluxa - 0.5*Fluxb - 0.5*Fluxc);
% Fluxbeta  = 2/3 * sqrt(3)/2 * (Fluxb - Fluxc);

% Fluxd = 2/3*( (Fluxa-0.5*Fluxb-0.5*Fluxc).*cos(theta) + sqrt(3)/2 * (Fluxb - Fluxc) .* sin(theta));
% Fluxq = 2/3*( -(Fluxa-0.5*Fluxb-0.5*Fluxc).*sin(theta) + sqrt(3)/2 * (Fluxb - Fluxc) .* cos(theta));

% Fluxd = Fluxd(2:(end));
% Fluxd = [Fluxd(NTime-1);Fluxd];
% Fluxq = Fluxq(2:(end));
% Fluxq = [Fluxq(NTime-1);Fluxq];
% 
% Fluxd = Fluxd(1:(end-1));
% Fluxq = Fluxq(1:(end-1));


% Fluxd_0 = mean(Fluxd);
% Fluxq_0 = mean(Fluxq);

Fluxa = zeros(size(theta,1),n3ph);
Fluxb = zeros(size(theta,1),n3ph);
Fluxc = zeros(size(theta,1),n3ph);
Fluxd = zeros(size(theta,1),n3ph);
Fluxq = zeros(size(theta,1),n3ph);

for ii=1:n3ph
    theta = theta+(th0(ii)-th0(1))*pi/180;
    
    Fluxa(:,ii) = Flux_1sim(:,1+3*(ii-1));
    Fluxb(:,ii) = Flux_1sim(:,2+3*(ii-1));
    Fluxc(:,ii) = Flux_1sim(:,3+3*(ii-1));
    
    Fluxd(:,ii) = 2/3*(+(Fluxa(:,ii)-0.5*Fluxb(:,ii)-0.5*Fluxc(:,ii)).*cos(theta)+sqrt(3)/2*(Fluxb(:,ii)-Fluxc(:,ii)).*sin(theta));
    Fluxq(:,ii) = 2/3*(-(Fluxa(:,ii)-0.5*Fluxb(:,ii)-0.5*Fluxc(:,ii)).*sin(theta)+sqrt(3)/2*(Fluxb(:,ii)-Fluxc(:,ii)).*cos(theta));
end

Fluxd = Fluxd(2:end,:);
Fluxd = [Fluxd(NTime-1,:);Fluxd];
Fluxq = Fluxq(2:end,:);
Fluxq = [Fluxq(NTime-1,:);Fluxq];

Fluxd = Fluxd(1:(end-1),:);
Fluxq = Fluxq(1:(end-1),:);

Fluxd = mean(Fluxd,2);
Fluxq = mean(Fluxq,2);





