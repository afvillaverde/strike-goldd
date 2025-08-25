clear all;

% States (1):
syms N
x = [N].';

% Model Parameters:
syms alpha K
p = [alpha K].';

% Known Inputs:
u = [];

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = [alpha*log(K/N)*N;
];

% Initial Conditions:
ics  = [];  

% Known Initial Conditions:
known_ics = [0];

% Outputs:
h = x;

% Saving:
save('GOM2','x','p','u','w','h','f','ics','known_ics');
