% Controllability example 30 from Vidyasagar.
clear;

% 2 states
syms x1 x2 x3
x = [x1; x2; x3];

% 1 input
syms u

% dynamic equations
f = [2*x1^2;
    x2^2 + x3*x2 + u*x3;
    -2*x3^3];

save('ctrl_ex_30_VIDYASAGAR','x','u','f');

% Controllability example 35 from Vidyasagar.
clear;

% 2 states
syms x1 x2 x3
x = [x1; x2; x3];

% 1 input
syms u

% dynamic equations
f = [u*(x1 + 2*x2 + 4*x3) - 14*x3;
    2*u*x2;
    3*u*x3 - 19*x3];

save('ctrl_ex_35_VIDYASAGAR','x','u','f');

% Controllability example 45 from Vidyasagar.
clear;

% 2 states
syms x1 x2 x3
x = [x1; x2; x3];

% 1 input
syms u

% dynamic equations
f = [3*x3 + u*(x1 + 2*x2 + 4*x3);
    6*x3 + u*(2*x1 + 2*x2);
    3*u*x3 - 2*x3];

save('ctrl_ex_45_VIDYASAGAR','x','u','f');