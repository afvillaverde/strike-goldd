% Markov model, analysed in:
% "Input-Dependent Structural Identifiability of Nonlinear Systems"
% AF Villaverde, ND Evans, MJ Chappell, JR Banga
% IEEE Control Systems Letters 3 (2), 272-277
% for analysis with u=constant and 4 experiments:
clear all;

% 2 -> 8 states
syms x1Exp1 x2Exp1 x1Exp2 x2Exp2 x1Exp3 x2Exp3 x1Exp4 x2Exp4
x = [x1Exp1 x2Exp1 x1Exp2 x2Exp2 x1Exp3 x2Exp3 x1Exp4 x2Exp4].';

% 1 -> 4 outputs
h = [x1Exp1 x1Exp2 x1Exp3 x1Exp4].';

% 1 -> 4 inputs
syms u1Exp1 u1Exp2 u1Exp3 u1Exp4
u = [u1Exp1 u1Exp2 u1Exp3 u1Exp4].';

% 10 unknown parameters 
syms a12 a21 b12 b21 a23 a32 b23 b32 a13 b13
p =[a12 a21 b12 b21 a23 a32 b23 b32 a13 b13].';

% dynamic equations
f = [x2Exp1*exp(a21 + b21*u1Exp1) - x1Exp1*(exp(a12 + b12*u1Exp1) + exp(a13 + b13*u1Exp1)) - exp(a13 - a12 + a21 - a23 + a32 + u1Exp1*(b13 - b12 + b21 - b23 + b32))*(x1Exp1 + x2Exp1 - 1); %Exp1
                                                   x1Exp1*exp(a12 + b12*u1Exp1) - exp(a32 + b32*u1Exp1)*(x1Exp1 + x2Exp1 - 1) - x2Exp1*(exp(a21 + b21*u1Exp1) + exp(a23 + b23*u1Exp1));     %Exp1
     x2Exp2*exp(a21 + b21*u1Exp2) - x1Exp2*(exp(a12 + b12*u1Exp2) + exp(a13 + b13*u1Exp2)) - exp(a13 - a12 + a21 - a23 + a32 + u1Exp2*(b13 - b12 + b21 - b23 + b32))*(x1Exp2 + x2Exp2 - 1); %Exp2
                                                   x1Exp2*exp(a12 + b12*u1Exp2) - exp(a32 + b32*u1Exp2)*(x1Exp2 + x2Exp2 - 1) - x2Exp2*(exp(a21 + b21*u1Exp2) + exp(a23 + b23*u1Exp2));     %Exp2
     x2Exp3*exp(a21 + b21*u1Exp3) - x1Exp3*(exp(a12 + b12*u1Exp3) + exp(a13 + b13*u1Exp3)) - exp(a13 - a12 + a21 - a23 + a32 + u1Exp3*(b13 - b12 + b21 - b23 + b32))*(x1Exp3 + x2Exp3 - 1); %Exp3
                                                   x1Exp3*exp(a12 + b12*u1Exp3) - exp(a32 + b32*u1Exp3)*(x1Exp3 + x2Exp3 - 1) - x2Exp3*(exp(a21 + b21*u1Exp3) + exp(a23 + b23*u1Exp3));     %Exp3
     x2Exp4*exp(a21 + b21*u1Exp4) - x1Exp4*(exp(a12 + b12*u1Exp4) + exp(a13 + b13*u1Exp4)) - exp(a13 - a12 + a21 - a23 + a32 + u1Exp4*(b13 - b12 + b21 - b23 + b32))*(x1Exp4 + x2Exp4 - 1); %Exp4
                                                   x1Exp4*exp(a12 + b12*u1Exp4) - exp(a32 + b32*u1Exp4)*(x1Exp4 + x2Exp4 - 1) - x2Exp4*(exp(a21 + b21*u1Exp4) + exp(a23 + b23*u1Exp4))];    %Exp4
% initial conditions
ics  = []; 
known_ics = [1,0,1,0,1,0,1,0];

save('Markov_3cycle_db_reduced_4Exp','x','p','h','f','u','ics','known_ics');