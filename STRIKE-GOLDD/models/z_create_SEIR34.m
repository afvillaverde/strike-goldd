clear;

% 5 states
syms N S E I R
x = [N; S; E; I; R];

% 7 unknown parameters
syms r beta mu epsilon gamma K
p = [r; beta; mu; epsilon; gamma; K];

% 2 outputs
h = [
    K*I;
    N
];

% one input
syms A
u = A;

% dynamic equations
f = [
    0;

    A - r*beta*S*I/N - mu*S;

    r*beta*S*I/N - epsilon*E - mu*E;

    epsilon*E - gamma*I - mu*I;

    gamma*I - mu*R
];

% initial conditions
ics = [];
known_ics = [];

save('SEIR34','x','p','h','f','u','ics','known_ics');