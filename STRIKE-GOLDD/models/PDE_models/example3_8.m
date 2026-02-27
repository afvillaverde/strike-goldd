%--------------------------------------------------------------------------
% File that creates the PDE model for fluorescence recovery after photobleaching 
% (FRAP). 
% The model is stored in a mat-file named example3_8.mat.
%--------------------------------------------------------------------------
% The model corresponds to Example presented in Section 5.1 of:
% Ciocanel MV, Ding L, Mastromatteo L, Reichheld S, Cabral S, Mowry K, 
% Sandstede B (2024) Parameter identifiability in PDE models of fluorescence 
% recovery after photobleaching.
%--------------------------------------------------------------------------

% PDE Model definition
clc; clear;

% Parameters and States(in terms of independent variables)
syms z(t,s) c(t,s) %states
syms D B1 B2 %parameters

% States and Observes states
x = {z c};
observed_vars = {z};
observation_eq = {};
% Parameters
p = {D B1 B2};

% PDE equations
z_t = diff(z,t);
c_t = diff(c,t);
c_ss = diff(c,s,2);
z_ss = diff(z,s,2);
z_t = D * z_ss - D * c_ss;
c_t = B2 * z - (B1+B2)* c;

f = {z_t;c_t};

%optional equations
opt_eq = {};

% save (Modelname.mat)
save('Example3_8.mat', 'x', 'observed_vars', 'p', 'f','observation_eq',"opt_eq");
