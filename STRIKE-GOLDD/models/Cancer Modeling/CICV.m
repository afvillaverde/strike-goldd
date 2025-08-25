clear all;

% States (2):
syms C I
x = [C; I];

% Model Parameters:
syms r1 r2 alpha1 k1 beta1 delta alpha2 alpha3 k2 beta2 gamma
p = [r1; r2; alpha1; k1; delta; alpha2; alpha3; k2; beta2];

% Known Inputs:
u = [beta1, gamma];

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = [
    r1*C*(1 - r2*C) - (alpha1*C*I)/(k1 + C) - beta1*C;
    delta - alpha2*C*I + (alpha3*C^2*I)/(C^2 + k2) - beta2*I + gamma
];

% Initial Conditions:
ics = [];

% Known Initial Conditions:
known_ics = [0, 0];

% Outputs:
h = [I];

% Saving:
save('CICV','x','p','u','w','h','f','ics','known_ics');
