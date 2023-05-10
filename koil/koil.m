% Data for the winding with Q = 12 slots and p = 5 pole pairs
% Winding factor:
kw = 9.330e-001;
% Coil pitch:   
yq = 1;

% slot matrix of phase A. Elements: 12
ka = [ 1.0, -0.5, 0.0, 0.0, 0.0, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0, -0.5]; 

% slot matrix of phase B. Elements: 12
kb = [ 0.0, 0.5, -1.0, 0.5, 0.0, 0.0, 0.0, -0.5, 1.0, -0.5, 0.0, 0.0]; 

% slot matrix of phase C. Elements: 12
kc = [ 0.0, 0.0, 0.0, -0.5, 1.0, -0.5, 0.0, 0.0, 0.0, 0.5, -1.0, 0.5]; 
