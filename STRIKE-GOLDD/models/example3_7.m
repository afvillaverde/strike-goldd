%--------------------------------------------------------------------------
% File that creates the general two-state reaction–diffusion PDE model. 
% It stores the model in a mat-file named example3_7.mat.
%--------------------------------------------------------------------------
% The model corresponds to Example (b) presented in Section 2 of:
% Browning AP, Taşca M, Falcó C, Baker RE (2024) Structural identifiability 
% analysis of linear reaction–advection–diffusion processes in mathematical 
% biology.
%--------------------------------------------------------------------------

% PDE Model definition
clc; clear;

% Parameters and States(in terms of independent variables)
syms q(t,s) v(t,s) %states
syms D1 D2 alpha1 alpha2 p1 p2 p3 p4 p5 p6 %parameters

% States and Observes states
x = {q v};
observed_vars = {};
observation_eq = {q+v};
% Parameters
p = {D1 D2 alpha1 alpha2 p1 p2 p3 p4 p5 p6};

% PDE equations
q_t = diff(q,t);
q_s = diff(q,s);
v_t = diff(v,t);
v_s = diff(v,s);
q_ss = diff(q,s,2);
v_ss = diff(v,s,2);
q_t = D1 * q_ss + alpha1 * q_s + p1 * q + p2 * v + p3;
v_t = D2 * v_ss + alpha2 * v_s + p4 * q + p5 * v + p6;

f = {q_t;v_t};

%optional equations
opt_eq = {};

% save (Modelname.mat)
save('Example3_7.mat', 'x', 'observed_vars', 'p', 'f','observation_eq',"opt_eq");
