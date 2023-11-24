% Model from J.W. Bolie. "Coefficients of normal blood glucose regulation". 
% J. Appl. Physiol., 16(5):783-788, 1961.

clear;

% 2 states
syms q1 q2
x = [q1; q2];

% 1 input
syms delta;
u = delta;

% 5 parameters 
syms p1 p2 p3 p4 Vp 
p =[p1; p2; p3; p4; Vp];

% 1 output
h = q1/Vp;

% initial conditions
ics  = []; 
known_ics = [0,0];

% dynamic equations
f = [-q1-2*q2+u;
    -3*q2+4*q1];

save('Bolie2','x','p','u','h','f','ics','known_ics');