%--------------------------------------------------------------------------
% Chemical reaction model.
% The model is taken from:
%--------------------------------------------------------------------------
% Benjamin Merkt et al. (2015) Higher-order Lie symmetries in 
% identifiability and predictability analysis of dynamic models
%--------------------------------------------------------------------------
clear all;

% 1 state:
syms x1 
x = [x1].';

% 3 parameters:
syms k s1 s2 
p = [k s1 s2].';

h = s1*x1/(1+s2*x1);

% No known inputs:
u = [];

% No unknown inputs:
w = [];

% dynamic equations:
f = [ 
	-2*x1*x1*k
];

% initial conditions:
ics  = [];   

% which initial conditions are known:
known_ics = [0]; 

save('CR','x','p','u','w','h','f','ics','known_ics');
