% 2-compartment linear model analysed in:
% "Input-Dependent Structural Identifiability of Nonlinear Systems"
% AF Villaverde, ND Evans, MJ Chappell, JR Banga
% IEEE Control Systems Letters 3 (2), 272-277
clear;

% 2 states
syms x1 x2
x = [x1; x2];

% 1 output
h = x1;

% one input
syms u1;
u = u1;

% 4 unknown parameters 
syms x3 x4 x5 x6
p =[x3; x4; x5; x6];

% dynamic equations
f = [-(x3+x4)*x1+x5*x2+x6*u1;
    x4*x1-x5*x2];

% initial conditions
ics  = []; 
known_ics = [1,0];

save('C2M','x','p','h','f','u','ics','known_ics');

