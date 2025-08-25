clear all;

% States (7):
syms C  C_d T_a A c4 p1 T_hat
x = [C; C_d; T_a; A; c4; p1; T_hat];

% Model Parameters:
syms lambda1 lambda2 pp phi a b rho psi sigma tau1 tau2 eta mu nu h s q iota
p = [lambda1 lambda2 pp phi a b rho psi sigma tau1 tau2 eta mu nu h s q iota].';

% Known Inputs:
syms K_C omega_i K_T delta_c4 delta_p1 i_c4 i_p1
u = [K_C; omega_i; K_T; delta_c4; delta_p1; i_c4; i_p1];

% Unknown inputs:
w = [];

% Dynamic Ecuations:
f = [
    lambda1*C*(1 - lambda2*(C + C_d)) - u(1) - pp*(1+p1)*((T_a/(C+C_d))^q / (s + (T_a/(C+C_d))^q))*C;
    u(1) - phi*u(2) - pp*(1+p1)*((T_a/(C+C_d))^q / (s + (T_a/(C+C_d))^q))*C_d;
    -u(3) + a*A - T_hat - iota*T_a*(C + C_d) - eta*T_a;
    rho*(C + C_d) + psi*phi*u(2)*C_d - sigma*A - a*A*T_hat - (b/(1+c4))*A*T_hat;
    u(6)*u(4) - nu*c4;
    u(7)*u(5) - mu*p1;
    -a*A*T_hat - (b/(1+c4))*A*T_hat + h
];

% Initial Conditions:
ics = [];

% Known Initial Conditions:
known_ics = [0, 0, 0, 0, 0, 0, 0];

% Outputs:
h = [C; A];

% Saving:
save('IMRAD','x','p','u','w','h','f','ics','known_ics');