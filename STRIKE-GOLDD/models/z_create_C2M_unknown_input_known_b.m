% 2-compartment linear model analysed in:
% "Input-Dependent Structural Identifiability of Nonlinear Systems"
% AF Villaverde, ND Evans, MJ Chappell, JR Banga
% IEEE Control Systems Letters 3 (2), 272-277
%--------------------------------------------------------------------------
% Modified so that 'b' is assumed to be known
 
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

% 4 unknown parameters 
syms x3 x4 x5 x6
p =[x3; x4; x5]; % ; x6 -> 'b' is a known constant

% dynamic equations
f = [-(x3+x4)*x1+x5*x2+x6*u1;
    x4*x1-x5*x2];

% initial conditions
ics  = []; 
known_ics = [0,0];

save('C2M_unknown_input_known_b','x','p','h','f','u','w','ics','known_ics');
