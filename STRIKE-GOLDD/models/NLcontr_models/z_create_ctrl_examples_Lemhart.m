% Controllability example 1 from Lenhart.
clear;

% 2 states
syms x1 x2
x = [x1; x2];

% 1 input
syms u

% dynamic equations
f = [u*x1;
    -x2*(u - 1)];

save('ctrl_ex_Lenhart_1','x','u','f');