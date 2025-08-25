clear all;

% States (2):
syms V_t N_t
x = [V_t; N_t];

% Model Parameters:
syms lambda K eta zeta gamma
p = [lambda K eta zeta gamma].';

% Known Inputs:
syms delta_ti
u = [delta_ti]; % Represents the combination of \delta(t-t_i)

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = [
    lambda * V_t * (1 - V_t/K) - eta * V_t - gamma * V_t * delta_ti;
    eta * V_t - zeta * N_t + gamma * V_t * delta_ti
];

% Initial Conditions:
ics = [];

% Known Initial Conditions:
known_ics = [0, 0];

% Outputs:
h = [V_t; N_t];

% Saving:
save('NECR','x','p','u','w','h','f','ics','known_ics');