% Sequential reaction model: A -> B -> C
clear;

% 3 states
syms xA xB xC
x = [xA; xB; xC];

% 2 outputs
h = [
    xC;
    xA + xB + xC
];

% y2(t) = eA * xA(t) + eB * xB(t) + eC * xC(t)

% no inputs
u = [];

% 2 unknown parameters
syms k1 k2
p = [k1; k2];

% dynamic equations
f = [
    -k1*xA;
     k1*xA - k2*xB;
     k2*xB
];

% initial conditions
ics = [];
known_ics = [];

save('Slowfast','x','p','h','f','u','ics','known_ics');