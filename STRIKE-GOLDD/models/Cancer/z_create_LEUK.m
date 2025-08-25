clear all;

% States (6):
syms T T_N E f_T h_T g_T
x = [T; T_N; E; f_T; h_T; g_T];

% Model Parameters:
syms p_T k_T p_NA p_AN p_E d_E g_MAX c_R c_F f_MIN m_K C_K
p = [p_T k_T p_NA p_AN p_E d_E g_MAX c_R c_F f_MIN m_K C_K].';

% Known Inputs:
syms e_t
u = [e_t];

% Unknown Inputs:
w = [];

% Dynamic Ecuations:
f = [
    p_T*T*(1 - T/k_T) - p_NA*T + p_AN*T_N - e_t*T - f_T*h_T*E;
    -p_AN*T_N + p_NA*T;
    p_E - d_E*E - g_T*E;
    (m_K*T)/(C_K + T) - f_T; % Dinámica de f_T
    (c_F^2 + f_MIN*T^2)/(c_F^2 + T^2) - h_T; % Dinámica de h_T
    (g_MAX*T^2)/(c_R^2 + T^2) - g_T % Dinámica de g_T
];

% Initial Conditions:
ics = [];

% Known Initial Conditions:
known_ics = [0, 0, 0, 0, 0, 0];

% Outputs:
h = [T; T_N; E; f_T; h_T; g_T];

% Saving:
save('LEUK','x','p','u','w','h','f','ics','known_ics');