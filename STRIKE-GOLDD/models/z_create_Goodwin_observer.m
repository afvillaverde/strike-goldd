% Goodwin model, as analysed in the high-gain observer paper.

clear;

% 3 states
syms x1 x2 x3 
x = [x1;x2;x3];

% 1 output
h = x1;

% no known input
u = [];

% 1 unknown input
syms f3
w = f3;

% known constants
syms b a AA k gamma delta

% unknown parameters
p = [];

% dynamic equations
f = [-b*x1 + a/(AA+k*x2);
    gamma*x3 - delta*x2;
    f3];

% initial conditions
ics  = []; 
known_ics = [0,0,0];

save('Goodwin_observer','x','p','h','u','w','f','ics','known_ics')
