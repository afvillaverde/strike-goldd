clear;

% 5 states
syms q1 q3 q35 q36 q7
x = [q1; q3; q35; q36; q7];

% one input
syms u
u = u;

% unknown parameters
syms k3 k4 k5 k6 k7 R V3 S V36
p = [k3; k4; k5; k6; k7; R; V3; S; V36];

% output
h = q7;

% dynamic equations
f = [
    k4*q3 - (k3 + k7)*q1 + u;

    k3*q1 - k4*q3 ...
    - k5*q3*(R*V3 - q35) + k6*q35 ...
    - k5*q3*(5*V36/V3)*(S*V36 - q36) + k6*q36;

    k5*q3*(R*V3 - q35) - k6*q35;

    k5*q3*(5*V36/V3)*(S*V36 - q36) - k6*q36;

    k7*q1
];

% initial conditions
ics = [];
known_ics = [];

save('CGV1990','x','p','h','f','u','ics','known_ics');