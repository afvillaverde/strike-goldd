% Controllability example 1.1 from Kawski (accessible but not controllable).
clear;

% 2 states
syms x1 x2
x = [x1; x2];

% 1 input
syms u

% no parameters 
p =[];

% output
h = x;

% initial conditions
ics  = []; 
known_ics = [0,0];

% dynamic equations
f = [u;
    x1^2];

save('ctrl_ex_11','x','p','u','h','f','ics','known_ics');

syms m
f = [u;
    x1^m];
save('ctrl_ex_31','x','p','u','h','f','ics','known_ics');
f = [u;
    x1^3];
save('ctrl_ex_31_odd','x','p','u','h','f','ics','known_ics');
f = [u;
    x1^4];
save('ctrl_ex_31_even','x','p','u','h','f','ics','known_ics');

syms x3
x = [x1; x2; x3];
f = [u;
    x1;
    x1^3*x2];
save('ctrl_ex_32','x','p','u','h','f','ics','known_ics');

f = [u;
    x1;
    x2^2+x1^3];
save('ctrl_ex_33','x','p','u','h','f','ics','known_ics');

syms x4
x = [x1; x2; x3; x4];
f = [u;
    x1;
    (1/6)*x1^3;
    x2*x3];
save('ctrl_ex_41','x','p','u','h','f','ics','known_ics');

syms lambda
f = [u;
    x1;
    x2;
    x2^2-lambda*x1^4];
save('ctrl_ex_51','x','p','u','h','f','ics','known_ics');

f = [u;
    x1;
    x2;
    x3^2-lambda*x1^4];
save('ctrl_ex_52','x','p','u','h','f','ics','known_ics');

f = [u;
    x1;
    x1^3;
    x3^2-lambda*x2^2*x1^4];
save('ctrl_ex_53','x','p','u','h','f','ics','known_ics');

f = [u;
    x1;
    x1^3;
    x3^2-x2^7];
save('ctrl_ex_61','x','p','u','h','f','ics','known_ics');

f = [u;
    x1;
    x1^3;
    x3^2-x2^8];
save('ctrl_ex_61b','x','p','u','h','f','ics','known_ics');


