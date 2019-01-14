% Model from J.W. Bolie. "Coefficients of normal blood glucose regulation". 
% J. Appl. Physiol., 16(5):783-788, 1961.
% modified so that it has only 3 unknown parameters (p1, p2 assumed known)

clear;

% 2 states
syms q1 q2
x = [q1; q2];

% 1 input
syms delta;
u = delta;

% 5 parameters, 3 unknown
syms p1 p2 p3 p4 Vp 
p =[p3; p4; Vp];

% 1 output
h = q1/Vp;

% initial conditions
ics  = []; 
known_ics = [0,0];

% dynamic equations
f = [p1*q1-p2*q2+u;
    p3*q2+p4*q1];

save('Bolie_C','x','p','u','h','f','ics','known_ics');