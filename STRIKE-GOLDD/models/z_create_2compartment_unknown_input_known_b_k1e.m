% 2-compartment linear model analysed in:
% "Input-Dependent Structural Identifiability of Nonlinear Systems"
% AF Villaverde, ND Evans, MJ Chappell, JR Banga
% IEEE Control Systems Letters 3 (2), 272-277
%--------------------------------------------------------------------------
% Modified so that 'b' and 'kie' are assumed to be known
clear;

% 2 states
syms x1 x2
x = [x1; x2];

% 1 output
h = x1;

% 0 known inputs
u = [];

% 1 unknown input
syms u1;
w = u1;

% 2 unknown parameters (x4 = k12, x5 = k21) 
syms x4 x5 
p =[x4; x5]; 

% 2 known constants (x3 = k1e, x6 = b)
syms x3 x6 

% dynamic equations
f = [-(x3+x4)*x1+x5*x2+x6*u1;
    x4*x1-x5*x2];

% initial conditions are in principle unknown
ics       = []; 
known_ics = [0,0];

save('two_compartment_unknown_input_known_b_k1e','x','p','h','f','u','w','ics','known_ics');
