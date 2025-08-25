clear all;

% RAD2 model adds K and theta.

% States (1):
syms V
x = [V];

% Model Parameters:
syms lambda alpha beta d K theta
p = [lambda alpha beta d K theta].';

% Known Inputs:
syms delta_ti
u = [delta_ti]; % Represents the combination of \delta(t-t_i)

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = lambda * V * (1 - (V/K)^theta) - (alpha*d + beta*d^2)*V*delta_ti;

% Initial Conditions:
ics = [];

% Known Initial Conditions:
known_ics = [0];

% Outputs:
h = [V];

% Saving:
save('RAD2','x','p','u','w','h','f','ics','known_ics');
