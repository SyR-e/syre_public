%%
% Assegnazione dei coefficienti dei materiali

% calore specifico carcassa (frame) [J/g/K]
Materiali.Carcassa.Cth=1005.6;
% Resistenza termica carcassa
Materiali.Carcassa.Rth=0.025;
% calore specifico flangia (Endcap) [J/kg/K]
Materiali.Carcassa.Cth=1005.6;
% densità di carcassa
Materiali.Carcassa.d=2700;
% Conducibilita del ferro di statore in laminazione assiale [W/m K] + MANUALE
Materiali.FerroSta.Rth=25;
% Conducibilità termica radiale del ferro di statore [W/m K]
Materiali.FerroSta.Rth=25;
% Densità del ferro
Materiali.FerroSta.d=7800;
% calore specifico del ferro di statore [J/kg/K]
Materiali.FerroSta.Cth=490;
Materiali.Conduttore.Rth=387;
% calore specifico del rame [J/kg/K]
Materiali.Conduttore.Cth=393;
% densità del rame [kg/m3]
Materiali.Conduttore.d=8500;
% coefficiente conduttività film di traferro (statico o dinamico) [W/m^2 K]
Materiali.Isolante.Rth=0.033;
% calore specifico aria
Materiali.Isolante.Cth=1014;
% conducibilità termica dell'acciaio dell'albero [W/m K]
Materiali.Albero.Rth=1;
% densità ferro albero
Materiali.Albero.d=7800;
% capacità termica del acciaio albero [J/kg/K]
Materiali.Albero.Cth=1;
% densità aria
Materiali.Isolante.d=1.2;
%
ManTherm.R1     = 0.0328; 
ManTherm.hc     = 1*10^3;

R1 = ManTherm.R1; 
hc = ManTherm.hc;

