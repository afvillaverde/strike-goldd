% 3-state linear model with one input
clear;

% 3 states
syms x1 x2 x3
x = [x1; x2; x3];

% 1 output
h = x1;

% one input
syms u0
u = u0;

% 5 unknown parameters
syms p1 p3 p4 p6 p7
p = [p1; p3; p4; p6; p7];

% dynamic equations
f = [
    -p1*x1 + x2 + u0;
     p3*x1 - p4*x2 + x3;
     p6*x1 - p7*x3
];

% initial conditions
ics = [];
known_ics = [];

save('Daisy_ex3','x','p','h','f','u','ics','known_ics');