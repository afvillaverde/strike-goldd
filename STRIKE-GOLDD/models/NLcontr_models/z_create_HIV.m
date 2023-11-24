%--------------------------------------------------------------------------
% HIV Model with Constant and Time Varying Parameters
% Miao H, Xia X, Perelson AS, Wu H. 
% "On identifiability of nonlinear ODE models and applications in viral dynamics." 
% SIAM review 53.1 (2011): 3-39.
%--------------------------------------------------------------------------
clear all;

% 3 states:
syms Tu Ti V
x = [Tu Ti V].';

% 5 parameters:
syms lambda rho N delta1 c 
p = [lambda rho N delta1 c].';
% 2 outputs:
h = [V, Tu+Ti].';

% 1 unknown input (time-varying parameter):
syms eta;
u = [];
w = eta;

% dynamic equations:
f = [lambda-rho*Tu-eta*Tu*V;
    eta*Tu*V-delta1*Ti;
    N*delta1*Ti-c*V];

% initial conditions:
ics  = [600,33,1e5];   

% which initial conditions are known:
known_ics = [0,0,1]; 

save('HIV','x','p','u','w','h','f','ics','known_ics');

u = eta;
w = [];
save('HIV_known_u','x','p','u','w','h','f','ics','known_ics');