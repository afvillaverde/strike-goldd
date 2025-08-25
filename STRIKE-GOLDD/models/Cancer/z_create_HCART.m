clear all;

% States (3):
syms CT CM T
x = [CT CM T].';

% Model Parameters:
syms phi rho epsilon theta alpha mu r b gamma
p = [phi rho epsilon theta alpha mu r b gamma].';

% Known Inputs:
u = [];

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = [
    phi * CT - rho * CT + theta * T * CM - alpha * T * CT;    % dCT/dt
    epsilon * CT - theta * T * CM - mu * CM;                  % dCM/dt
    r * T * (1 - b * T) - gamma * CT * T;
];

% Initial Conditions:
ics  = [];  

% Known Initial Conditions:
known_ics = [0, 0, 0];

% Outputs:
h = [T];

% Saving:
save('HCART','x','p','u','w','h','f','ics','known_ics');
