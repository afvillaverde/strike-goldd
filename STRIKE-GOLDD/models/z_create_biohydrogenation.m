% Model with states x4, x5, x6, x7
clear;

% 4 states
syms x4 x5 x6 x7
x = [x4; x5; x6; x7];

% 2 outputs
h = [x4; x5];

% no inputs
u = [];

% 6 unknown parameters
syms k5 k6 k7 k8 k9 k10
p = [k5; k6; k7; k8; k9; k10];

% dynamic equations
f = [
    -k5*x4/(k6 + x4);
    
     k5*x4/(k6 + x4) ...
     - k7*x5/(k8 + x5 + x6);
    
     k7*x5/(k8 + x5 + x6) ...
     - k9*x6*(k10 - x6)/k10;
    
     k9*x6*(k10 - x6)/k10
];

% initial conditions
ics = [];
known_ics = [];

save('Biohydrogenation','x','p','h','f','u','ics','known_ics');