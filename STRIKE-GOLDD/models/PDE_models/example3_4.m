%--------------------------------------------------------------------------
% File that creates the reaction–diffusion PDE model of cancer invasion. 
% It stores the model in a mat-file named example3_4.mat.
%--------------------------------------------------------------------------
% The model corresponds to Example 15 presented in Section 5 of:
% Byrne HM, Harrington HA, Ovchinnikov A, Pogudin G, Rahkooy H, Soto P 
% (2025) Algebraic identifiability of partial differential equation models.
%--------------------------------------------------------------------------

% PDE Model definition
clc; clear;

% Parameters and States(in terms of independent variables)
syms q(t,s) v(t,s) w(t,s) %states
syms r1 r2 r3 k1 k2 d1 d2 d3 d4 a12 a21  %parameters

% States and Observes states
x = {q v w};
observed_vars = {q v w};
observation_eq = {};
% Parameters
p = {r1 r2 r3 k1 k2 d1 d2 d3 d4 a12 a21};

% PDE equations
q_t = diff(q,t);
v_t = diff(v,t);
w_t = diff(w,t);
v_s = diff(v,s);

w_ss = diff(w,s,2);

q_t = r1*q*(1- q/k1 -a12*v/k2) -d1*w*q;
v_t = r2*v*(1- v/k2 - a21*q/k1)+ d2*diff(((1-q)*v_s),s);
w_t = d4*w_ss + r3*v - d3*w;
f = {q_t;v_t;w_t};

%optional equations
opt_eq = {};

% save (Modelname.mat)
save('Example3_4.mat', 'x', 'observed_vars', 'p', 'f','observation_eq','opt_eq');
