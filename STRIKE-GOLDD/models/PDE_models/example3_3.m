%--------------------------------------------------------------------------
% File that creates the coupled reaction–diffusion PDE model. 
% It stores the model in a mat-file named example3_3.mat.
%--------------------------------------------------------------------------
% The model corresponds to Example 13 presented in Section 5 of:
% Byrne HM, Harrington HA, Ovchinnikov A, Pogudin G, Rahkooy H, Soto P 
% (2025) Algebraic identifiability of partial differential equation models.
%--------------------------------------------------------------------------

% PDE Model definition
clc; clear;

% Parameters and States(in terms of independent variables)
syms q(t,s) v(t,s) %states
syms d1 d2 a1 a2 b1 b2 c1 c2 %parameters

% States and Observes states
x = {q v};
observed_vars = {q v};
observation_eq = {};
% Parameters
p = {d1 d2 a1 a2 b1 b2 c1 c2};

% PDE equations
q_t = diff(q,t);
v_t = diff(v,t);
q_ss = diff(q,s,2);
v_ss = diff(v,s,2);
q_t = d1*q_ss + q*(a1 - b1*q - c1*v);
v_t = d2*v_ss + v*(a2 - b2*q - c2*v);

f = {q_t;v_t};

%optional equations
opt_eq = {};

% save (Modelname.mat)
save('Example3_3.mat', 'x', 'observed_vars', 'p', 'f','observation_eq','opt_eq');
