%--------------------------------------------------------------------------
% File that creates the linear reaction–advection–diffusion PDE model 
% describing a two-state cell cycle model of cell migration subject to 
% exponential growth. The model is stored in a mat-file named example3_6.mat.
%--------------------------------------------------------------------------
% The model corresponds to Example (a) presented in Section 2 of:
% Browning AP, Taşca M, Falcó C, Baker RE (2024) Structural identifiability 
% analysis of linear reaction–advection–diffusion processes in mathematical 
% biology.
%--------------------------------------------------------------------------

% PDE Model definition
clc; clear;

% Parameters and States(in terms of independent variables)
syms r(t,s) g(t,s)
syms D1 D2 k1 k2

% States and Observes states
x = {r g};
observed_vars = {};
observation_eq = {r+g};
% Parameters
p = {D1 D2 k1 k2};

% PDE equations
g_ss = diff(g,s,2);
r_ss = diff(r,s,2);

r_t = D1 * r_ss - k1*r + 2*k2*g;
g_t = D2 * g_ss + k1*r - k2*g;

f = {r_t; g_t};

%optional equations
opt_eq = {};

% save (Modelname.mat)
save('Example3_6.mat', 'x', 'observed_vars', 'p', 'f','observation_eq',"opt_eq");
