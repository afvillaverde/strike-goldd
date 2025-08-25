clear all;

% States (8):
syms C_I C_E C_P T_P T_N M_i M_a IL_6
x = [C_I C_E C_P T_P T_N M_i M_a IL_6].';

% Model Parameters:
syms eta mu_I nu kappa epsilon theta mu_E mu_P rho K gama A B C...
   sigma_M delta_M sigma_I alpha delta_I g0 beta_B beta_K beta_C
p = [eta mu_I nu kappa epsilon theta mu_E mu_P rho K gama A B C...
   sigma_M delta_M sigma_I alpha delta_I g0 beta_B beta_K beta_C].';

% Known Inputs:
u = [];

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = [
  - eta * (T_P / (A + T_P)) * C_I - mu_I * C_I;
  nu * (T_P / (A + T_P)) * C_I + kappa * (T_P / (A + T_P)) * C_E - epsilon * (1 - (T_P / (A + T_P))) * C_E + theta * (T_P / (A + T_P)) * C_P - mu_E * C_E;
  epsilon * (1 - (T_P / (A + T_P))) * C_E - theta * (T_P / (A + T_P)) * C_P - mu_P * C_P;
  rho * T_P * (1 - (T_P + T_N) / K) - gama * (C_E / (B + C_E)) * T_P;
  rho * T_N * (1 - (T_P + T_N) / K) - g0 * gama * (C_E / (B + C_E)) * T_N;
  sigma_M - (beta_B*(T_P/(A+T_P))*C_E+beta_K*(C_E/(B+C_E))*(T_P+g0*T_N)+beta_C*(M_a/(C+M_a))*C_E)* M_i - delta_M * M_i;
  (beta_B*(T_P/(A+T_P))*C_E+beta_K*(C_E/(B+C_E))*(T_P+g0*T_N)+beta_C*(M_a/(C+M_a))*C_E)* M_i - delta_M * M_a;
  sigma_I + alpha * M_a - delta_I * IL_6;
];

% Initial Conditions:
ics = [];

% Known Initial Conditions:
known_ics = [0, 0, 0, 0, 0, 0, 0, 0];

% Outputs:    
h = [C_E, M_i + M_a, IL_6];

% Saving:
save('CRS','x','p','u','w','h','f','ics','known_ics');
