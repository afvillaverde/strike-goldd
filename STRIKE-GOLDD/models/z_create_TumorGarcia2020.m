% Tumor-effector model
clear;

% 2 states
syms T E
x = [T; E];

% 1 output
h = E;

% one input
syms sigma
u = sigma;

% 5 unknown parameters
syms a b k d m
p = [a; b; k; d; m];

% dynamic equations
f = [
    a*T*(1 - b*T) - k*T*E;      % Tumor cells
    sigma - d*E + m*E*T         % Effector cells
];

% initial conditions
ics = [];
known_ics = [];

save('TumorGarcia2020','x','p','h','f','u','ics','known_ics');