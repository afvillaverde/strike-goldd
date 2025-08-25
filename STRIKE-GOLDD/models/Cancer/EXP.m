clear all;

% States (1):
syms V
x = [V].';

% Model Parameters:
syms lambda K theta
p = [lambda, K, theta].';

% Known Inputs:
u = [];

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = [(lambda/theta)*V*(1-(V/K)^theta);
];

% Initial Conditions:
ics  = [];  

% Known Initial Conditions:
known_ics = [0];

% Output:
h = x;

% Saving:
save('EXP','x','p','u','w','h','f','ics','known_ics');
