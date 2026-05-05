% 3-compartment linear model with one input
clear;

% 3 states
syms x1 x2 x3
x = [x1; x2; x3];

% 1 output
h = x1;

% one input
syms u
u = u;

% 6 unknown parameters
syms a21 a31 a01 a12 a13
p = [a21; a31; a01; a12; a13];

% dynamic equations
f = [
    -(a21 + a31 + a01)*x1 + a12*x2 + a13*x3 + u;
     a21*x1 - a12*x2;
     a31*x1 - a13*x3
];

% initial conditions
ics = [];
known_ics = [];

save('Daisy_Mamil3','x','p','h','f','u','ics','known_ics');