clear all;

% States (2):
syms V_S V_R
x = [V_S V_R].';

% Model Parameters:
syms lambda_S lambda_R K_S K_R gamma_R gamma_S  
p = [lambda_S, lambda_R, K_S, K_R, gamma_R, gamma_S].';

% Known Inputs:
u = [];

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = [lambda_S*V_S*(1-(V_S/K_S)-gamma_R*(V_R/K_S));
    lambda_R*V_R*(1-(V_R/K_R)-gamma_S*(V_S/K_R));
];

% Initial Conditions:
ics  = [];  

% Known Initial Conditions:
known_ics = [0, 0];

% Outputs:
h = V_S+V_R;

% Saving:
save('L_V','x','p','u','w','h','f','ics','known_ics');
