clear all;

% States (1):
syms V
x = [V].';

% Model Parameters:
syms lambda
p = [lambda].';

% Known Inputs:
u = [];

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = [lambda*V;
];

% Initial Conditions:
ics  = [];  

% Known Initial Conditions:
known_ics = [0];

% Output:
h = x;

% Saving:
save('EXP1','x','p','u','w','h','f','ics','known_ics');
