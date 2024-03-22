% Data for the winding with Q = 12 slots and p = 2 pole pairs
% Winding factor:
kw = 1.000e+000;
% Coil pitch:   
yq = 3;

% slot matrix of phase A. Elements: 12
ka = [ 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0]; 

% slot matrix of phase B. Elements: 12
kb = [ 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0]; 

% slot matrix of phase C. Elements: 12
kc = [ 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0]; 
