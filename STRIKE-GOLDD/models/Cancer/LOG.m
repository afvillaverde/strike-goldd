clear all;

% States (1):
syms N
x = [N].';

% Model Parameters:
syms a K
p = [a K].';

% Known Inputs:
u = [];

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = [a*N*(1-(N/K));
];

% Initial Conditions:
ics  = [];  

% Known Initial Conditions:
known_ics = [0];

% Outputs:
h = x;

% Saving:
save('LOG','x','p','u','w','h','f','ics','known_ics');
