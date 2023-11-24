% CSTR model used in the high-gain observer paper.

clear;

% 3 states
syms Cr Tr Tc 
x = [Cr;Tr;Tc];

% 1 output
h = Tc;

% 3 known inputs:
syms u1 u2 u3
u = [u1;u2;u3];

% 1 unknown input
syms d
w = d;

% known constants
syms q qc gamma

% unknown parameters
syms delta1
p = delta1;

% aux functions
F = exp(gamma*Tr/(gamma+Tr/20));

% dynamic equations
f = [q*(u1-Cr)-0.072*F*Cr+d;
     q*(u2-Tr)+0.576*F*Cr-0.3*(Tr-Tc);
     delta1*qc*(u3-Tc)+3.0*(Tr-Tc)];

% initial conditions
ics  = []; 
known_ics = [0,0,0];

save('CSTR_observer','x','p','h','u','w','f','ics','known_ics')
