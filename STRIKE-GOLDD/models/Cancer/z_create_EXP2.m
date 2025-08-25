clear all;

% States (1):
syms V
x = [V].';

% Model Parameters:
syms lambda K
p = [lambda, K].';

% Known Inputs:
u = [];

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = [lambda*V*(1-(V/K));
];

% Initial Conditions:
ics  = [];  

% Known Initial Conditions:
known_ics = [0];

% Output:
h = x;

% Saving:
save('EXP2','x','p','u','w','h','f','ics','known_ics');
