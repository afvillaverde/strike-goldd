% SEIR model with constant population N
clear;

% 5 states
syms S E I R
x = [S; E; I; R];

% no inputs
u = [];

% 3 unknown parameters
syms beta alpha lambda N
p = [beta; alpha; lambda];

% 2 outputs
h = [
    I
];

% dynamic equations
f = [
    -beta*I*(S/N);

     beta*I*(S/N) - alpha*E;

     alpha*E - lambda*I;

     lambda*I;
];

% initial conditions
ics = [];
known_ics = [];

save('SEIRT','x','p','h','f','u','ics','known_ics');