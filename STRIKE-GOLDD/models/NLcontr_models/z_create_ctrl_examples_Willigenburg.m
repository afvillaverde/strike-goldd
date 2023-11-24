% Controllability example 1 from Willigenburg.
clear;

% 2 states
syms x1 x2 x3 x4
x = [x1; x2; x3; x4];

% 1 input
syms u1 u2
u = [u1;u2];

% dynamic equations
f = [u1*cos(x3 + x4);
    u1*sin(x3 + x4);
    u1*sin(x4);
    u2];

save('ctrl_ex_1_Willigenburg','x','u','f');

% Controllability example 2 from Willigenburg.
clear;

% 2 states
syms x1 x2 x3
x = [x1; x2; x3];

% 1 input
syms u1 u2
u = [u1;u2];

% dynamic equations
f = [u1*x2^2;
    u2*x3^2;
    x1^2 - 1];

save('ctrl_ex_2_Willigenburg','x','u','f');