%% VERSIONE 20 11 2011
% Voltage Elab

% dq e.m.f.
Ed = per.EvalSpeed * geo.p * pi/30 * (-Fluxq);
Eq = per.EvalSpeed * geo.p * pi/30 * ( Fluxd);

% Ed_g = per.EvalSpeed * geo.p * pi/30 * ( Fluxq);
% Eq_g = per.EvalSpeed * geo.p * pi/30 * ( Fluxd);

pos_EMF = (pos_Mov(2:end) + pos_Mov(1:end-1))/2;
dT = (time(2) - time(1))/1000;
% Ea_ =  (Fluxa(2:end)-Fluxa(1:end-1))/dT;
% Eb_ =  (Fluxb(2:end)-Fluxb(1:end-1))/dT;
% Ec_ =  (Fluxc(2:end)-Fluxc(1:end-1))/dT;

% Ea = interp1(pos_EMF,Ea_,pos_Mov,'pchip');
% Eb = interp1(pos_EMF,Eb_,pos_Mov,'pchip');
% Ec = interp1(pos_EMF,Ec_,pos_Mov,'pchip');

% Ea = Ea(1:(end-1),:);
% Eb = Eb(1:(end-1),:);
% Ec = Ec(1:(end-1),:);


Ea_ =  (Fluxa(2:end,:)-Fluxa(1:end-1,:))/dT;
Eb_ =  (Fluxb(2:end,:)-Fluxb(1:end-1,:))/dT;
Ec_ =  (Fluxc(2:end,:)-Fluxc(1:end-1,:))/dT;

Ea = zeros(length(time),n3ph);
Eb = zeros(length(time),n3ph);
Ec = zeros(length(time),n3ph);

for ii=1:n3ph
    Ea(:,ii) = interp1(pos_EMF,Ea_(:,1),pos_Mov,'pchip');
    Eb(:,ii) = interp1(pos_EMF,Eb_(:,1),pos_Mov,'pchip');
    Ec(:,ii) = interp1(pos_EMF,Ec_(:,1),pos_Mov,'pchip');
end

Ea = Ea(1:(end-1),:);
Eb = Eb(1:(end-1),:);
Ec = Ec(1:(end-1),:);
