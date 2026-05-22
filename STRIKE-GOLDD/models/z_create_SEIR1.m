% SEIQR model
clear;

% 5 states
syms S E I R Q
x = [S; E; I; R; Q];

% no inputs
u = [];

% unknown parameters
syms beta v psi gamma
p = [beta; v; psi; gamma];

% outputs
h = [
    Q
    ];

% dynamic equations
f = [
    beta*S*I;
    beta*S*I - v*E;
    v*E - psi*I - (1 - psi)*gamma*I;
    gamma*Q + (1 - psi)*gamma*I;
    -gamma*Q + psi*I;
    ];

% initial conditions
ics = [];
known_ics = [];

save('SEIR1','x','p','h','f','u','ics','known_ics');