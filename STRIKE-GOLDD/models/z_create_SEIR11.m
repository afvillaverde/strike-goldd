clear;

% 6 states
syms S E I R F
x = [S; E; I; R; F];

% no inputs
u = [];

% 5 unknown parameters
syms beta xi sigma gamma mu N
p = [beta; xi; sigma; gamma; mu];

% 2 outputs
h = [
    F
];

% dynamic equations
f = [

    -beta*S*I/N + xi*R;

     beta*S*I/N - sigma*E;

     sigma*E - gamma*I - mu*I;

     gamma*I - xi*R;

     mu*I
];

% initial conditions
ics = [];
known_ics = [];

save('SEIR11','x','p','h','f','u','ics','known_ics');