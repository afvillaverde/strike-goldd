%--------------------------------------------------------------------------
% File that creates the model of the age-structured SEI epidemic PDE model. 
% It stores it in a mat-file named SEI.mat.
%--------------------------------------------------------------------------
% The model corresponds to the example presented in Section 4.2 of:
% Renardy M, Kirschner D, Eisenberg M (2022) Structural identifiability 
% analysis of age-structured PDE epidemic models. 
%--------------------------------------------------------------------------

% PDE Model definition
clc; clear;
% Parameters and States(in terms of independent variables)
syms S(t,s) yE(t,s) yI(t,s) %states
syms beta kI muS epsilon kE delta muE muI %parameters

% States and Observes states
x = {S yE yI};
observed_vars = {yE yI};
observation_eq = {};
% Parameters
p = {beta kI muS epsilon kE delta muE muI};

% PDE equations
S_t = diff(S,t);
S_s = diff(S,s);
yE_t = diff(yE,t);
yE_s = diff(yE,s);
yI_t = diff(yI,t);
yI_s = diff(yI,s);

S_t = -S_s -beta*S*yI/kI - muS*S;
yE_t = -yE_s + (1-epsilon)*beta*s*yI*kE/kI - (delta+muE)* yE;
yI_t = -yI_s + epsilon*beta*s*yI + delta*yE*kI/kE - muI*yI ;

f = {S_t; yE_t; yI_t};

%optional equations
syms u_muS u_muE
opt_eq = {u_muS - u_muE ==0};

% save (Modelname.mat)
save('Example3_1.mat', 'x', 'observed_vars', 'p', 'f','observation_eq','opt_eq');
