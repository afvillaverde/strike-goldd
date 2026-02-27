%--------------------------------------------------------------------------
% File that creates the PDE model corresponding to the well-known wave 
% equation. 
% It stores the model in a mat-file named example3_5.mat.
%--------------------------------------------------------------------------
% The model corresponds to Example 17 presented in Section 5 of:
% Byrne HM, Harrington HA, Ovchinnikov A, Pogudin G, Rahkooy H, Soto P 
% (2025) Algebraic identifiability of partial differential equation models. 
%--------------------------------------------------------------------------

% PDE Model definition
clc; clear;

% Parameters and States(in terms of independent variables)
syms q(t,s)  %states
syms delta  %parameters

% States and Observes states
x = {q };
observed_vars = {q};
observation_eq = {};
% Parameters
p = {delta};

% PDE equations
q_tt = diff(q,t,2);
q_ss = diff(q,s,2);
q_tt = delta^2 * q_ss;

f = {q_tt};

%optional equations
opt_eq = {};

% save (Modelname.mat)
save('Example3_5.mat', 'x', 'observed_vars', 'p', 'f','observation_eq',"opt_eq");
