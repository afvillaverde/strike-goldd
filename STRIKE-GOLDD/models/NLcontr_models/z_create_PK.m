%--------------------------------------------------------------------------
% Pharmacokinetic (PK) model.
% The model is taken from:
%--------------------------------------------------------------------------
% Merkt et al (2015) "Higher-order Lie symmetries in identifiability and
% predictability analysis of dynamic models" Phys Rev E 92(1-1):012920.
%--------------------------------------------------------------------------
clear all;

% 4 states:
syms x1 x2 x3 x4
x = [x1 x2 x3 x4].';

% 10 parameters:
syms k1 k2 k3 k4 k5 k6 k7...
    s2 s3
p = [k1 k2 k3 k4 k5 k6 k7...
    s2 s3].';

% 2 outputs:
h = [s2*x2; s3*x3];

% 1 input:
syms u1 ;

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

u = u1;
w = [];
save('PK_known_input','x','p','u','w','h','f','ics','known_ics');
u = [];
w = u1;
save('PK_unknown_input','x','p','u','w','h','f','ics','known_ics');