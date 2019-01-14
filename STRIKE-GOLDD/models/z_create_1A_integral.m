% Hormonal circuit model with integral feedback
% Originally published in: Karin et al, Mol Syst Biol 2016
% Corresponds to the model in Fig. 1A of: Villaverde & Banga, arXiv:1701.02562

clear;

% 2 states
syms x1  x2
x = [x1; x2];

% 1 output
h = x1;

% 1 input
syms u u0;

% 2 parameters 
syms p1 p2    
p =[p1; p2];

% initial conditions
syms x10 x20
ics  = [x10 x20]; 
known_ics = [0,0];

% dynamic equations
f = [u0+u-p2*x1-p1*x2;
    x1-x10];

save('1A_integral','x','p','h','u','f','ics','known_ics');