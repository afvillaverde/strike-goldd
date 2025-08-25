clear all;

% States (1):
syms N
x = [N].';

% Model Parameters:
syms a b gamma
p = [a b gamma].';

% Known Inputs:
u = [];

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = [a*N^gamma-b*N;
];

% Initial Conditions:
ics  = [];  

% Known Initial Conditions:
known_ics = [0];

% Outputs:
h = x;

% Saving:
save('BERT','x','p','u','w','h','f','ics','known_ics');
