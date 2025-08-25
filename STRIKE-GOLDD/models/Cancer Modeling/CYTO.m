clear all;

% This version adds phi_dot and ddphi, as external known inputs.

% States (2):
syms V dV
x = [V; dV];

% Model Parameters:
syms a
p = [a];

% Known Inputs:
syms phi_dot ddphi b
u = [phi_dot; ddphi; b];

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = [
    dV;
    (2/V)*dV^2 - (phi_dot*(a - b) - ddphi/phi_dot)*dV - (phi_dot)^2 * a * b * V
];

% Initial Conditions:
ics = [];

% Known Initial Conditions:
known_ics = [0, 0];

% Outputs:
h = [V];

% Saving:
save('CYTO2','x','p','u','w','h','f','ics','known_ics');

