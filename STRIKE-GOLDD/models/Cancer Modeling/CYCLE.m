clear all;

% States (4):
syms x1 x2 n V
x = [x1; x2; n; V];

% Model Parameters:
syms AG AR a epsilon theta delta alphaG K
p = [AG; AR; a; epsilon; theta; delta; alphaG; K];

% Known Inputs:
syms C(t) theta(t)
u = [C(t); theta(t)];

% Unknown Inputs:
w = [];

% Dynamic equations:
f = [
    (AG * x2 - AR * x1) * x1;
    (C(t) + (1 - C(t)) * AR * x2) * x2 - (AG * x2 - AR * x1) * x2;
    epsilon * n * (1 - n) * (theta(t) * x2 + (x2 - 1));
    (delta + alphaG) * V * (1 - V / K);
];

% Initial Conditions:
ics = [];

% Known Initial Conditions:
known_ics = [0, 0, 0, 0];

% Outputs:
h = V; 

% Saving:
save('CYCLE','x','p','u','h','f','ics','known_ics');
