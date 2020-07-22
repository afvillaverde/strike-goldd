%--------------------------------------------------------------------------
% Pharmacokinetic (PK) model.
% The model is taken from:
%--------------------------------------------------------------------------
% Benjamin Merkt et al. (2015) Higher-order Lie symmetries in 
% identifiability and predictability analysis of dynamic models
%--------------------------------------------------------------------------
clear all;

% 4 states:
syms x1 x2 x3 x4
x = [x1 x2 x3 x4].';

% 10 parameters:
syms k1 k2 k3 k4 k5 k6 k7...
    s2 s3...
    x2ob x3ob 
p = [k2 k3 k4 k5 k6 k7...
    k1 s2 s3].';

% 2 outputs:
x2ob=s2*x2;
x3ob=s3*x3;
h = [ x2ob x3ob ].';

% 1 known input:
syms u1 ;
u = u1;

% 0 unknown inputs:
w = [];

% dynamic equations:
f = [ 
	u1-(k1+k2)*x1;
    k1*x1-(k3+k6+k7)*x2+k5*x4;
    k2*x1+k3*x2-k4*x3;
    k6*x2-k5*x4
];

% initial conditions:
ics  = [];   

% which initial conditions are known:
known_ics = [0,0,0,0]; 

save('PK','x','p','u','w','h','f','ics','known_ics');