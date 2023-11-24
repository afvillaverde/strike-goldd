% Controllability example from Sussmann 1983.
clear;

% 2 states
syms x1 x2 x3
x = [x1; x2; x3];

% 1 input
syms u

% dynamic equations
f = [u;
    x1;
    x1^3 + x2^2];

save('ctrl_ex_Sussmann83','x','u','f');

% Controllability example from Sussmann 1987.
clear;

% 2 states
syms x1 x2 x3
x = [x1; x2;x3];

% 1 input
syms u

% dynamic equations
f = [u;
    x1;
    x1^3*x2];

save('ctrl_ex_Sussmann87','x','u','f');