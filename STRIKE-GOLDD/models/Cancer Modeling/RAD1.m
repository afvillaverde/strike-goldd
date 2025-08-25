clear all;

% States (1):
syms V
x = [V];

% Model Parameters:
syms lambda alpha beta d
p = [lambda alpha beta d].';

% Known Inputs:
syms delta_ti
u = [delta_ti]; % Represents the combination of \delta(t-t_i)

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = lambda * V - (alpha*d + beta*d^2)*V*delta_ti;

% Initial Conditions:
ics = [];

% Known Initial Conditions:
known_ics = [0];

% Outputs:
h = [V];

% Saving:
save('RAD1','x','p','u','w','h','f','ics','known_ics');