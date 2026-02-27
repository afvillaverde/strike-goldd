%--------------------------------------------------------------------------
% File that creates the PDE model. 
% It stores the model in a mat-file named Example3_2.mat
%--------------------------------------------------------------------------
% The model is taken from Example 7 in Section 3 of:
% Byrne HM, Harrington HA, Ovchinnikov A, Pogudin G, Rahkooy H, Soto P 
% (2025) Algebraic identifiability of partial differential equation models. 
%--------------------------------------------------------------------------

% PDE Model definition
clc; clear;

% Parameters and States(in terms of independent variables)
syms q(t,s) v(t,s) %states
syms a b c %parameters

% States and Observes states
x = {q v};
observed_vars = {v};
observation_eq = {};
% Parameters
p = {a b c};

% PDE equations
q_t = diff(q,t);
q_s = diff(q,s);
v_t = diff(v,t);
q_ss = diff(q,s,2);
q_t = a * q_s + c * v^2;
v_t = b * q_ss;

f = {q_t;v_t};

%optional equations
opt_eq = {}; 

% save (Modelname.mat)
save('Example3_2.mat', 'x', 'observed_vars', 'p', 'f','observation_eq','opt_eq');
