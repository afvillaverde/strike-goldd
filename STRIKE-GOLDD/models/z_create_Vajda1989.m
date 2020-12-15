% Unidentifiable nonlinear model
% S. Vajda et al. (1989) Similarity transformation approach to 
% identifiability analysis of nonlinear compartmental models

clear;

% 2 states
syms x1 x2 
x = [x1; x2];

% 1 unknown input
syms w 
u = [];

% 5 parameters 
syms t1 t2 t3 t4 
p =[t1;t2;t3;t4];

% 1 output
h = x1;

% initial conditions
ics  = [0,0]; 
known_ics = [1,1];

% dynamic equations
f = [t1*x1^2+t2*x1*x2+w;
    t3*x1^2+t4*x1*x2];

save('Vajda1989','x','p','w','u','h','f','ics','known_ics');