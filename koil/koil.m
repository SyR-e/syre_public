% Data for the winding with Q = 24 slots and p = 2 pole pairs
% Winding factor:
kw = 9.659e-001;
% Coil pitch:   
yq = 6;

% slot matrix of phase A. Elements: 24
ka = [ 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0]; 

% slot matrix of phase B. Elements: 24
kb = [ 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0]; 

% slot matrix of phase C. Elements: 24
kc = [ 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]; 
