function [Nut,NuEw] = NusseltNumb(DinamicoStatico)
% calcolo del numero di Nusselt

global PM PC

%% calcolo dei numeri per il traferro
% calcolo dei numeri di Taylor e Prandtl TODO 
Nta = 1;
Npr = 1;

if DinamicoStatico    
    if (Nta <= 41)
        Nut = 2.2;
    else
        Nut = 0.23 * (Nta^0.63) * (Npr^0.27);
    end
else
    Nut = 2.0;
end

%% calcolo del numero di Nusselt per le connessionni frontali
% numero di Reynolds e di Rayleigh TODO
Nre = 1;
Nra = 1;

NuEw=0.83*(Nre*Nra)^0.2;
