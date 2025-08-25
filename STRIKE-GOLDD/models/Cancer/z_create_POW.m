clear all;

% States (1):
syms N
x = [N].';

% Model Parameters:
syms a gamma
p = [a, gamma].';

% Known Inputs:
u = [];

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = [a*N^gamma;
];

% Initial Conditions:
ics  = [];  

% Known Initial Conditions:
known_ics = [0];

% Outputs:
h = x;

% Saving:
save('POW','x','p','u','w','h','f','ics','known_ics');
